#ifndef CIE_GL_UTILITY_GL_ERROR_HPP
#define CIE_GL_UTILITY_GL_ERROR_HPP

// --- Utility Includes ---
#include "packages/logging/inc/Loggee.hpp"

// --- STL Includes ---
#include <string>


namespace cie::gl {


void checkGLErrors( const std::string& r_message = "GL Error" );

void checkGLErrors( utils::Loggee& r_loggee,
                    const std::string& r_message = "GL Error" );


} // namespace cie::gl

#include "packages/utility/impl/GLError_impl.hpp"

#endif