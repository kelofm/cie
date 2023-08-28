#ifndef CIE_GL_TRAITS_TYPES_HPP
#define CIE_GL_TRAITS_TYPES_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- STL Includes ---
#include <string>


namespace cie::gl {


GLenum getGLTypeEnum( const std::string& r_typeName );

template <class T>
GLenum getGLTypeEnum();


} // namespace cie::gl

#include "packages/traits/impl/types_impl.hpp"

#endif