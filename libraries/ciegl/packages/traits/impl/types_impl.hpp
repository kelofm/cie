#ifndef CIE_GL_TRAITS_TYPES_IMPL_HPP
#define CIE_GL_TRAITS_TYPES_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <class T>
inline GLenum
getGLTypeEnum()
{
    CIE_THROW( Exception, "Unrecognized GL type" )
}


template <>
inline GLenum
getGLTypeEnum<char>()
{
    return GL_BYTE;
}


template <>
inline GLenum
getGLTypeEnum<unsigned char>()
{
    return GL_UNSIGNED_BYTE;
}


template <>
inline GLenum
getGLTypeEnum<short>()
{
    return GL_UNSIGNED_SHORT;
}


template <>
inline GLenum
getGLTypeEnum<unsigned short>()
{
    return GL_UNSIGNED_SHORT;
}


template <>
inline GLenum
getGLTypeEnum<int>()
{
    return GL_INT;
}


template <>
inline GLenum
getGLTypeEnum<unsigned int>()
{
    return GL_UNSIGNED_INT;
}


//template <>
//inline GLenum
//getGLTypeEnum<half>()
//{
//    return GL_HALF_FLOAT;
//}


template <>
inline GLenum
getGLTypeEnum<float>()
{
    return GL_FLOAT;
}


template <>
inline GLenum
getGLTypeEnum<double>()
{
    return GL_DOUBLE;
}


} // namespace cie::gl


#endif