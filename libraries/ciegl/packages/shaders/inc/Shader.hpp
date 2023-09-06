#ifndef CIE_GL_SHADER_HPP
#define CIE_GL_SHADER_HPP

// --- External Includes ---
#include <glad/glad.h>

// --- Utility Includes ---
#include "packages/types/inc/IDObject.hpp"
#include "packages/types/inc/NamedObject.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <string>
#include <iostream>
#include <memory>
#include <filesystem>


namespace cie::gl {


 /** @brief Interface class representing a shader written in GLSL.
  *
  *  @details This is the common base class for vertex, geometry, fragment,
  *           and compute shaders. Construction requires the GLSL source
  *           as well as a configuration in JSON format, defining the
  *           various in- and outputs of the shader. The configuration is
  *           expected to conform to the following format:
  *           {
  *               "attributes" : {
  *                   "<attribute name>" {
  *                       "size"      : <number of components in the attribute>,
  *                       "stride"    : <number of components between two attribute instances>,
  *                       "offset"    : <number of components before the first component of the first attribute>,
  *                       "type"      : "<data type>"
  *                   },
  *                   ...
  *               },
  *               "uniforms" : {
  *                   "<uniform name>" : {
  *                       "size"  : <number of components in the uniform>,
  *                       "type"  : "<data type>"
  *                   },
  *                   ...
  *               },
  *               "textures" : {
  *                   "<texture name>" : {
  *                       "dimensions" : <number of texture dimensions [1, 2, 3] (eg.: 2 for images)>,
  *                       "channels" : <number of components each entry in the texture consists of (eg.: 3 for RGB images)>,
  *                       "type" : <data type of the texture's components>
  *                   },
  *                   ...
  *               },
  *               "outputs" : {
  *                   "TODO" : ""
  *               }
  *           }
  *  @ingroup ciegl
  */
class Shader : public utils::IDObject<Size>
{
public:
    class Attribute;

    class Uniform;

    class Texture;

    class Output;

    using AttributeContainer = DynamicArray<Attribute>;

    using UniformContainer = DynamicArray<Uniform>;

    using TextureContainer = DynamicArray<Texture>;

    using OutputContainer = DynamicArray<Output>;

    CIE_DEFINE_CLASS_POINTERS(Shader)

public:
    CIE_DEFINE_CLASS_DEFAULT_MOVES(Shader)

    Shader(const Shader&) = delete;

    Shader& operator=(const Shader&) = delete;

    /// @brief Destructor calling @p glDeleteShader.
    virtual ~Shader();

    /// @brief Print the GLSL source to the specified stream.
    /// @param r_stream: output stream to print the source to.
    friend Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const Shader> r_shader);

    /// @brief Load the GLSL source on demand.
    virtual std::string source() const = 0;

    /// @brief Load the JSON configuration on demand.
    virtual std::string configuration() const = 0;

    /// @brief Get the set of registered attributes parsed from the JSON configuration file.
    const AttributeContainer& attributes() const;

    /// @brief Get the set of registered uniforms parsed from the JSON configuration file.
    const UniformContainer& uniforms() const;

    /// @brief Get the set of registered textures parsed from the JSON configuration file.
    const TextureContainer& textures() const;

    /// @brief Get the set of registered outputs parsed from the JSON configuration file.
    const OutputContainer& outputs() const;

public:
    /** @brief Immutable attribute descriptor for vertex shaders.
     *  @details Represents an attribute passed to each vertex.
     *           For example: position, texture coordinates, etc.
     *           Requires:
     *            - name   : name of the variable in the shader (eg: position)
     *            - size   : number of components that make up the attribute (eg: x,y,z -> 3)
     *            - stride : number of bytes in between the components of the attribute
     *                       (that are not part of the attribute itself)
     *            - offset : number of bytes before the first component of the attribute
     */
    class Attribute : public utils::NamedObject
    {
    public:
        /** @brief Attribute constructor.
         *
         *  @param r_name: attribute name in the GLSL source.
         *  @param size: number of components that make up the attribute (eg.: x, y, z => size 3)
         *  @param stride: number of bytes in between the components of the attribute
         *                 (that are not part of the attribute itself)
         *  @param offset: number of bytes before the first component of the first attribute.
         */
        Attribute(const std::string& r_name,
                  Size size,
                  Size stride,
                  Size offset) :
            utils::NamedObject(r_name),
            _size(size),
            _stride(stride),
            _offset(offset)
        {}

        /// Size getter.
        Size size() const { return _size; }

        /// Stride getter.
        Size stride() const { return _stride; }

        /// Offset getter.
        Size offset() const { return _offset; }

    private:
        const Size _size;

        const Size _stride;

        const Size _offset;
    };


    /** @brief Immutable uniform descriptor (for any shader).
     *  @details Represents data global to all entities of the shader type
     *           (vertices for vertex shaders, fragments for fragment shaders, etc.).
     *           For example: transformation matrix, background color, etc.
     *           Requires:
     *            - name : name of the uniform in the shader code
     *            - size : number of components that make up the uniform
     *            - type : data type
     */
    class Uniform : public utils::NamedObject
    {
    public:
        /** @brief Uniform Constructor.
         *
         *  @param r_name: name of the uniform in the GLSL source.
         *  @param size: number of components that make up the uniform (eg.: x, y, t => size 3).
         *  @param r_type: type of the uniform (eg.: vec3, mat4, etc.).
         *  @todo: fix ugly uniform binding.
         */
        Uniform(const std::string& r_name,
                 Size size,
                 const std::string& r_type) :
            utils::NamedObject(r_name),
            _size(size),
            _type(r_type)
        {}

        Size size() const         { return _size; }

        const std::string& type() const { return _type; }

    protected:
        const Size _size;

        const std::string _type;
    };


    /** @brief Immutable texture descriptor.
     *  @details Represents a texture in a GLSL shader.
     *           Requires:
     *            - name: name of the texture in the GLSL source.
     *            - dimension: texture dimension (1, 2, or 3).
     *            - channels: number of components per pixel in the texture (eg.: RGB => 3).
     *            - type: data type of the texture's components.
     */
    class Texture : public utils::NamedObject
    {
    public:
        /** @brief Texture Constructor.
         *  @param r_name name of the texture in the GLSL source.
         *  @param dimension texture dimension (1, 2, or 3).
         *  @param channels number of components per pixel in the texture (eg.: RGB => 3).
         *  @param r_type data type of the texture's components.
         */
        Texture(const std::string& r_name,
                 Size dimension,
                 Size channels,
                 const std::string& r_type) :
            utils::NamedObject(r_name),
            _dimension(dimension),
            _channels(channels),
            _type(r_type)
        {}

        Size dimension() const {return _dimension;}

        Size channels() const {return _channels;}

        const std::string& type() const {return _type;}

    protected:
        Size _dimension;

        Size _channels;

        std::string _type;
    };

    /// Represents an output of the fragment shader.
    class Output :
        public utils::NamedObject,
        public utils::IDObject<Size>
    {
    public:
        Output(const std::string& r_name, Size id) :
            utils::NamedObject(r_name),
            utils::IDObject<Size>(id)
        {}
    };

protected:
    Shader() noexcept;

    void parseConfig(Ref<const std::string> r_configuration);

private:
    AttributeContainer _attributes;

    UniformContainer _uniforms;

    TextureContainer _textures;

    OutputContainer _outputs;
}; // class Shader



/// @brief @ref Shader loaded from the filesystem.
class DynamicShader final : public Shader
{
public:
    /// @copydoc Shader::source
    std::string source() const override;

    /// @copydoc Shader::configuration
    std::string configuration() const override;

    /// @brief Get the GLSL source file's path.
    Ref<const std::filesystem::path> sourcePath() const;

    /// @brief Get the JSON configuration file's path.
    Ref<const std::filesystem::path> configurationPath() const;

private:
    /** @brief Construct a shader by directly specifying its source and configuration files.
     *  @param r_configPath: path to the JSON configuration file.
     *  @param r_sourcePath: path to the GLSL source file.
     */
    DynamicShader(const std::filesystem::path& r_configPath,
                  const std::filesystem::path& r_sourcePath);

    /** @brief Construct a shader by name and containing directory.
     *  @details two files are expected:
     *            - <@p r_directory>/@p r_name.glsl (shader)
     *            - <@p r_directory>/@p r_name.json (configuration)
     *  @param r_name: name of the shader, or matching shader set.
     *  @param r_directory: path to the directory containing the GLSL and the JSON configuration.
     */
    DynamicShader(const std::string& r_name,
                  const std::filesystem::path& r_directory);

    friend Shader::SharedPointer shaderFactory(const std::filesystem::path& r_configPath,
                                               const std::filesystem::path& r_sourcePath,
                                               GLenum shaderType);

private:
    std::filesystem::path _sourcePath;

    std::filesystem::path _configPath;
}; // class DynamicShader



/// @brief @ref Shader loaded from memory.
class StaticShader final : public Shader
{
public:
    std::string source() const override;

    std::string configuration() const override;

private:
    StaticShader(Ref<const std::string> r_key);

    friend Shader::SharedPointer shaderFactory(Ref<const std::string> r_key,
                                               GLenum shaderType);

private:
    Ptr<const std::string> _p_source;

    Ptr<const std::string> _p_configuration;
}; // class StaticShader


/** @brief Vertex shader factory.
 *  @details Create a vertex shader from a GLSL source file and its
 *           JSON configuration file. The configuration is checked for
 *           the basic requirements of a vertex shader.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeVertexShader(const std::filesystem::path& r_configPath,
                                       const std::filesystem::path& r_codePath);


/** @brief Vertex shader factory.
 *  @details Create a vertex shader from its GLSL source and json configuration
 *           loaded from the library.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeVertexShader(Ref<const std::string> r_key);


/** @brief Geometry shader factory.
 *  @details Create a geometry shader from a GLSL source file and its
 *           JSON configuration file. The configuration is checked for
 *           the basic requirements of a geometry shader.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeGeometryShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_codePath);


/** @brief Geometry shader factory.
 *  @details Create a geometry shader from its GLSL source and json configuration
 *           loaded from the library.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeGeometryShader(Ref<const std::string> r_key);


/** @brief Fragment shader factory.
 *  @details Create a fragment shader from a GLSL source file and its
 *           JSON configuration file. The configuration is checked for
 *           the basic requirements of a fragment shader.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeFragmentShader(const std::filesystem::path& r_configPath,
                                         const std::filesystem::path& r_codePath);


/** @brief Fragment shader factory.
 *  @details Create a fragment shader from its GLSL source and json configuration
 *           loaded from the library.
 *  @note The consistency between the source and configuration is not checked.
 *  @ingroup ciegl
 */
Shader::SharedPointer makeFragmentShader(Ref<const std::string> r_key);


} // namespace cie::gl


#endif
