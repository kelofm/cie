#ifndef CIE_GL_UTILITY_GL_ERROR_IMPL_HPP
#define CIE_GL_UTILITY_GL_ERROR_IMPL_HPP

// --- External Includes ---
#include <glad/glad.h>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


inline void checkGLErrors( const std::string& r_message )
{
    GLint id = glGetError();

    if ( id ) [[unlikely]]
        CIE_THROW(
            Exception,
            r_message + " | error ID: " + std::to_string(id)
        )
}


inline void checkGLErrors( utils::Loggee& r_loggee,
                           const std::string& r_message )
{
    GLint id = glGetError();

    if ( id ) [[unlikely]]
        r_loggee.log( r_message + " (error " + std::to_string(id) + ")",
                      LOG_TYPE_ERROR );
}


} // namespace cie::gl


#endif