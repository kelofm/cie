// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"
#include <packages/types/inc/IDObject.hpp>

// --- GL Includes ---
#include "packages/shaders/inc/Shader.hpp"
#include "ciegl_shaders.hpp"

// --- STL Includes ---
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>


namespace cie::gl {


Shader::Shader() noexcept
    : utils::IDObject<Size>(std::numeric_limits<Size>::max())
{
}


void Shader::parseConfig(Ref<const std::string> r_configuration)
{
    CIE_BEGIN_EXCEPTION_TRACING
    io::JSONObject configuration(r_configuration);

    // Parse config - attributes
    {
        const auto attributes = configuration["attributes"];
        const auto it_end = attributes.end();
        for (auto it=attributes.begin(); it!=it_end; ++it)
        {
            const auto item = *it;
            _attributes.emplace_back(
                it.key(),
                item["size"].as<Size>(),
                item["stride"].as<Size>(),
                item["offset"].as<Size>()
            );
        }
    }

    // Parse config - uniforms
    {
        const auto uniforms = configuration["uniforms"];
        const auto it_end = uniforms.end();
        for (auto it=uniforms.begin(); it!=it_end; ++it)
        {
            const auto item = *it;
            _uniforms.emplace_back(
                it.key(),
                item["size"].as<Size>(),
                item["type"].as<std::string>()
            );
        }
    }

    // Parse config - textures
    {
        const auto textures = configuration["textures"];
        const auto it_end = textures.end();
        for (auto it=textures.begin(); it!=it_end; ++it)
        {
            const auto item = *it;
            _textures.emplace_back(
                it.key(),
                item["dimension"].as<Size>(),
                item["channels"].as<Size>(),
                item["type"].as<std::string>()
            );
        }
    }

    // Parse config - outputs
    {
        const auto outputs = configuration["outputs"];
        const auto it_end = outputs.end();
        Size index = 0;
        for (auto it=outputs.begin(); it!=it_end; ++it)
        {
            _outputs.emplace_back(
                it.key(),
                index++
            );
        }
    }
    CIE_END_EXCEPTION_TRACING
}


DynamicShader::DynamicShader(const std::filesystem::path& r_configPath,
                             const std::filesystem::path& r_sourcePath)
    : Shader(),
      _sourcePath(r_sourcePath),
      _configPath(r_configPath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (std::filesystem::is_symlink(_sourcePath))
        _sourcePath = std::filesystem::read_symlink(_sourcePath);

    if (std::filesystem::is_symlink(_configPath))
        _configPath = std::filesystem::read_symlink(_configPath);

    parseConfig(this->configuration());

    CIE_END_EXCEPTION_TRACING
}



Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const Shader> r_shader)
{
    r_stream << r_shader.configuration() << '\n' << r_shader.source() << '\n';
    return r_stream;
}


Shader::~Shader()
{
    glDeleteShader(this->getID());
}


const typename Shader::AttributeContainer& Shader::attributes() const
{
    return _attributes;
}


const typename Shader::UniformContainer& Shader::uniforms() const
{
    return _uniforms;
}


const typename Shader::TextureContainer& Shader::textures() const
{
    return _textures;
}


const typename Shader::OutputContainer& Shader::outputs() const
{
    return _outputs;
}


std::string DynamicShader::source() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(std::filesystem::is_regular_file(_sourcePath),
              "File '" << _sourcePath << "' is not a file")

    std::ifstream source(_sourcePath);
    CIE_CHECK(source.is_open(), "Unable to open shader source file: '" << _sourcePath.string() << "'")

    std::stringstream stream;
    stream << source.rdbuf();
    return stream.str();

    CIE_END_EXCEPTION_TRACING
}


std::string DynamicShader::configuration() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(std::filesystem::is_regular_file(_configPath),
              "File '" << _configPath << "' is not a file")

    std::ifstream source(_configPath);
    CIE_CHECK(source.is_open(), "Unable to open shader configuration file: '" << _configPath.string() << "'")

    std::stringstream stream;
    stream << source.rdbuf();
    return stream.str();

    CIE_END_EXCEPTION_TRACING
}


Ref<const std::filesystem::path> DynamicShader::sourcePath() const
{
    return _sourcePath;
}


Ref<const std::filesystem::path> DynamicShader::configurationPath() const
{
    return _configPath;
}


StaticShader::StaticShader(Ref<const std::string> r_key)
    : Shader()
{
    const auto& definition = ShaderDefinitions::get(r_key);
    _p_source = &definition.source;
    _p_configuration = &definition.configuration;
    this->parseConfig(this->configuration());
}


std::string StaticShader::source() const
{
    return *_p_source;
}


std::string StaticShader::configuration() const
{
    return *_p_configuration;
}


namespace impl {


void compileShader(Ref<Shader> r_shader, GLenum shaderType)
{
    CIE_BEGIN_EXCEPTION_TRACING
    // Check shader type
    // TODO: update when adding support for compute shaders
    CIE_CHECK(shaderType == GL_VERTEX_SHADER || shaderType == GL_GEOMETRY_SHADER || shaderType == GL_FRAGMENT_SHADER,
              "Unsupported shader type: " << shaderType)

    const GLuint id = glCreateShader(shaderType);
    r_shader.setID(Size(id));

    CIE_CHECK(glIsShader(r_shader.getID()) == GL_TRUE,
              "Failed to create " << (shaderType==GL_VERTEX_SHADER ? "vertex" : (shaderType==GL_GEOMETRY_SHADER) ? "geometry" : (shaderType==GL_FRAGMENT_SHADER) ? "fragment" : "unknown") << " shader with ID " << r_shader.getID())

    // Feed the source to OpenGL
    const auto source = r_shader.source();
    const char* p_source = source.c_str();
    glShaderSource(id, 1, &p_source, NULL);

    // Compile and check for errors
    glCompileShader(id);
    GLint compileStatus;
    glGetShaderiv(id, GL_COMPILE_STATUS, &compileStatus);

    if (compileStatus != GL_TRUE) [[unlikely]] {
        GLchar buffer[512];
        glGetShaderInfoLog(id, sizeof(buffer)/sizeof(GLchar), NULL, buffer);

        std::stringstream message;
        message << "Shader '" << id
        << "' failed to compile with message: " << std::endl
        << buffer << std::endl << source;
        CIE_THROW(Exception, message.str())
    }
    CIE_END_EXCEPTION_TRACING
}


void checkVertexShader(Ref<const Shader> r_shader)
{
    // Vertex shader needs at least one attribute
    CIE_CHECK(!r_shader.attributes().empty(),
              "Vertex shader " << r_shader.getID() << " has no attributes\n"
              << r_shader.configuration() << '\n'
              << r_shader.source())
}


void checkGeometryShader(Ref<const Shader> r_shader)
{
}


void checkFragmentShader(Ref<const Shader> r_shader)
{
    // Vertex shader needs at least one attribute
    CIE_CHECK(!r_shader.outputs().empty(),
              "Fragment shader " << r_shader.getID() << " has no outputs\n"
              << r_shader.configuration() << '\n'
              << r_shader.source())
}


} // namespace impl


Shader::SharedPointer shaderFactory(const std::filesystem::path& r_configPath,
                                    const std::filesystem::path& r_sourcePath,
                                    GLenum shaderType)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = Shader::SharedPointer(
        new DynamicShader(r_configPath, r_sourcePath)
    );
    impl::compileShader(*p_shader, shaderType);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer shaderFactory(Ref<const std::string> r_key, GLenum shaderType)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = Shader::SharedPointer(
        new StaticShader(r_key)
    );
    impl::compileShader(*p_shader, shaderType);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeVertexShader(const std::filesystem::path& r_configPath,
                                       const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_configPath, r_sourcePath, GL_VERTEX_SHADER);
    impl::checkVertexShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeVertexShader(Ref<const std::string> r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_key, GL_VERTEX_SHADER);
    impl::checkVertexShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeGeometryShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_configPath, r_sourcePath, GL_GEOMETRY_SHADER);
    impl::checkGeometryShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeGeometryShader(Ref<const std::string> r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_key, GL_GEOMETRY_SHADER);
    impl::checkGeometryShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeFragmentShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_configPath, r_sourcePath, GL_FRAGMENT_SHADER);
    impl::checkFragmentShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeFragmentShader(Ref<const std::string> r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING
    auto p_shader = shaderFactory(r_key, GL_FRAGMENT_SHADER);
    impl::checkFragmentShader(*p_shader);
    return p_shader;
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl
