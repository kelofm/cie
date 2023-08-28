#ifndef CIE_GL_TRAITS_MATRIX_TYPES_IMPL_HPP
#define CIE_GL_TRAITS_MATRIX_TYPES_IMPL_HPP

// --- External Includes ---
#include "glm/glm.hpp"


namespace cie::gl {


template <>
inline GLenum
getGLTypeEnum<glm::mat2>()
{
    return GL_FLOAT_MAT2;
}


template <>
inline GLenum
getGLTypeEnum<glm::mat3>()
{
    return GL_FLOAT_MAT3;
}


template <>
inline GLenum
getGLTypeEnum<glm::mat4>()
{
    return GL_FLOAT_MAT4;
}


} // namespace cie::gl


#endif