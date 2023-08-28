// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"

// --- GL Includes ---
#include "packages/shaders/inc/Shader.hpp"

// --- STL Includes ---
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>


namespace cie::gl {


Shader::Shader(const std::filesystem::path& r_configPath,
               const std::filesystem::path& r_sourcePath) :
    utils::IDObject<Size>(std::numeric_limits<Size>().max()),
    _sourcePath(r_sourcePath),
    _configPath(r_configPath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    {
        if (std::filesystem::is_symlink(_sourcePath))
            _sourcePath = std::filesystem::read_symlink(_sourcePath);

        if (std::filesystem::is_symlink(_configPath))
            _configPath = std::filesystem::read_symlink(_configPath);

        io::JSONObject configuration(_configPath);

        // Parse config file - attributes
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

        // Parse config file - uniforms
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

        // Parse config file - textures
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

        // Parse config file - outputs
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
    }

    CIE_END_EXCEPTION_TRACING
}


Shader::~Shader()
{
    glDeleteShader(this->getID());
}


const std::filesystem::path& Shader::sourcePath() const
{
    return _sourcePath;
}


const std::filesystem::path& Shader::configurationPath() const
{
    return _configPath;
}


std::string Shader::source() const
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


std::string Shader::configuration() const
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


const typename Shader::attribute_container& Shader::attributes() const
{
    return _attributes;
}


const typename Shader::uniform_container& Shader::uniforms() const
{
    return _uniforms;
}


const typename Shader::texture_container& Shader::textures() const
{
    return _textures;
}


const typename Shader::output_container& Shader::outputs() const
{
    return _outputs;
}


Shader::SharedPointer shaderFactory(const std::filesystem::path& r_configPath,
                                    const std::filesystem::path& r_sourcePath,
                                    GLenum shaderType)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Construct shader base
    auto p_shader = Shader::SharedPointer(
        new Shader(r_configPath, r_sourcePath)
    );

    // Check shader type
    // TODO: update when adding support for compute shaders
    CIE_CHECK(shaderType == GL_VERTEX_SHADER || shaderType == GL_GEOMETRY_SHADER || shaderType == GL_FRAGMENT_SHADER,
              "Shader loaded from '" << r_sourcePath << "' and configured from '" << r_configPath
              << "' given unsupported type '" << shaderType << "'")

    const GLuint id = glCreateShader(shaderType);
    p_shader->setID(Size(id));

    CIE_CHECK(glIsShader(p_shader->getID()) == GL_TRUE,
              "Failed to create " << (shaderType==GL_VERTEX_SHADER ? "vertex" : (shaderType==GL_GEOMETRY_SHADER) ? "geometry" : (shaderType==GL_FRAGMENT_SHADER) ? "fragment" : "unknown") << " shader with ID " << p_shader->getID())

    // Feed the source to OpenGL
    const auto source = p_shader->source();
    const char* p_source = source.c_str();
    glShaderSource(id, 1, &p_source, NULL);

    // Compile and check for errors
    glCompileShader(id);
    GLint compileStatus;
    glGetShaderiv(id, GL_COMPILE_STATUS, &compileStatus);

    if (compileStatus != GL_TRUE) [[unlikely]]
    {
        GLchar buffer[512];
        glGetShaderInfoLog(id, sizeof(buffer)/sizeof(GLchar), NULL, buffer);

        std::stringstream message;
        message << "Shader '" << id << "' loaded from '" << r_sourcePath
        << "' and configured from '" << r_configPath
        << "' failed to compile with message: " << std::endl
        << buffer << std::endl << source;
        CIE_THROW(Exception, message.str())
    }

    return p_shader;

    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeVertexShader(const std::filesystem::path& r_configPath,
                                       const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_shader = shaderFactory(r_configPath, r_sourcePath, GL_VERTEX_SHADER);

    // Vertex shader needs at least one attribute
    CIE_CHECK(!p_shader->attributes().empty(),
              "Vertex shader " << p_shader->getID() << " loaded from '" << p_shader->sourcePath()
              << "' and configured from '" << p_shader->configurationPath() << "' has no attributes")

    return p_shader;

    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeGeometryShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Geometry shaders have no configuration requirements.
    return shaderFactory(r_configPath, r_sourcePath, GL_GEOMETRY_SHADER);

    CIE_END_EXCEPTION_TRACING
}


Shader::SharedPointer makeFragmentShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_sourcePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_shader = shaderFactory(r_configPath, r_sourcePath, GL_FRAGMENT_SHADER);

    // Vertex shader needs at least one attribute
    CIE_CHECK(!p_shader->outputs().empty(),
              "Vertex shader " << p_shader->getID() << " loaded from '" << p_shader->sourcePath()
              << "' and configured from '" << p_shader->configurationPath() << "' has no outputs")

    return p_shader;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl