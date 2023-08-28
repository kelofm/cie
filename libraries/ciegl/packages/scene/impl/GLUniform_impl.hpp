#ifndef CIE_GL_GL_UNIFORM_IMPL_HPP
#define CIE_GL_GL_UNIFORM_IMPL_HPP

// --- External Includes ---
#include <glm/gtc/type_ptr.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <class BoundType>
FloatVec3GLUniform<BoundType>::FloatVec3GLUniform( const Shader::Uniform& r_shaderUniform,
                                                   Size programID,
                                                   const BoundType& r_uniform ) :
    GLUniform( r_shaderUniform, programID ),
    _r_uniform( r_uniform )
{
    CIE_CHECK(
        this->_properties.type == GL_FLOAT_VEC3,
        "Uniform type mismatch: " + this->name()
        + " | type: " + std::to_string(this->_properties.type)
    )
}


template <class BoundType>
inline void
FloatVec3GLUniform<BoundType>::write()
{
    CIE_BEGIN_EXCEPTION_TRACING

    typename FloatVec3GLUniform<BoundType>::uniform_type uniform = _r_uniform;
    glUniform3fv( this->getID(),
                  1,
                  glm::value_ptr(uniform) );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif