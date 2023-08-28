#ifndef CIE_GL_TRAITS_VECTOR_TYPES_IMPL_HPP
#define CIE_GL_TRAITS_VECTOR_TYPES_IMPL_HPP

// --- External Includes ---
#include "glm/glm.hpp"


namespace cie::gl {


template <>
inline GLenum
getGLTypeEnum<glm::vec2>()
{
    return GL_FLOAT_VEC2;
}


template <>
inline GLenum
getGLTypeEnum<glm::vec3>()
{
    return GL_FLOAT_VEC3;
}


template <>
inline GLenum
getGLTypeEnum<glm::vec4>()
{
    return GL_FLOAT_VEC4;
}


} // namespace cie::gl


#endif