#ifndef CIE_GL_GL_UNIFORM_HPP
#define CIE_GL_GL_UNIFORM_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <glm/glm.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/IDObject.hpp"

// --- Internal Includes ---
#include "packages/shaders/inc/Shader.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::gl {




/* --- GLUniform --- */

/// Interface for loading uniforms into GL.
class GLUniform :
    public Shader::Uniform,
    public utils::IDObject<GLint>
{
private:
    struct Properties
    {
        GLsizei length;
        GLsizei size;
        GLenum  type;
    };

public:
    CIE_DEFINE_CLASS_POINTERS(GLUniform)

public:
    /** @brief Construct a uniform from a shader's uniform data, specific to a program.
     *
     *  @param r_shaderUniform uniform data loaded from a shader
     *  @param programID id of the program containing the shader
     */
    GLUniform(const Shader::Uniform& r_shaderUniform, Size programID);

    /// @brief Upload the uniform's data to the GPU.
    virtual void write() = 0;

    const Properties& properties() const;

protected:
    Properties _properties;
};


using GLUniformPtr = std::shared_ptr<GLUniform>;




/* --- GLUniformPlaceholder --- */

/**
 * Class for storing uniform information from shaders
 * during the cunstruction of a scene. After construction,
 * instances of this class must be replaced by specialized
 * GLUniforms of the correct type by binding specific objects
 * to them.
 */
class GLUniformPlaceholder final : public GLUniform
{
public:
    GLUniformPlaceholder( const Shader::Uniform& r_shaderUniform,
                          Size programID );
    void write() override;
};




/* --- FloatMat4GLUniform --- */

/**
 * 4x4 float matrix uniform used mostly for 3D transformations.
 * Should be bound to the transformation matrix of a camera.
 */
class FloatMat4GLUniform final : public GLUniform
{
public:
    using uniform_type = glm::mat4;

public:
    FloatMat4GLUniform( const Shader::Uniform& r_shaderUniform,
                        Size programID,
                        const uniform_type& r_uniform );

    void write() override;

private:
    const uniform_type& _r_uniform;
};



/* --- FloatVec3GLUniform --- */

/**
 * 3-component float vector uniform used mostly for the position of
 * cameras, light sources, or some sort of general direction.
 *
 * The bound object is copied at every updating to allow conversion
 * between different component types
 */
template <class BoundType>
class FloatVec3GLUniform final : public GLUniform
{
public:
    using bound_type    = BoundType;
    using uniform_type  = glm::vec3;

public:
    FloatVec3GLUniform( const Shader::Uniform& r_shaderUniform,
                        Size programID,
                        const bound_type& r_uniform );

    void write() override;

private:
    const bound_type& _r_uniform;
};



/* --- FloatGLUniform --- */

/**
 * 4x4 float matrix uniform used mostly for 3D transformations.
 * Should be bound to the transformation matrix of a camera.
 */
class FloatGLUniform final : public GLUniform
{
public:
    using uniform_type = GLfloat;

public:
    FloatGLUniform( const Shader::Uniform& r_shaderUniform,
                    Size programID,
                    const uniform_type& r_uniform );

    void write() override;

private:
    const uniform_type& _r_uniform;
};


} // namespace cie::gl

#include "packages/scene/impl/GLUniform_impl.hpp"

#endif