// --- Internal Includes ---
#include "packages/context/inc/GLFWContext.hpp"

// --- STL Includes ---
#include <iostream>
#include <string>

namespace cie::gl {


void GLAPIENTRY messageCallback(    GLenum source,
                                    GLenum type,
                                    GLuint id,
                                    GLenum severity,
                                    GLsizei length,
                                    const GLchar* message,
                                    const void* logger )
{
    std::string messageString = message;
    auto p_context = const_cast<GLFWContext*>( reinterpret_cast<const GLFWContext*>( logger ) );
    p_context->logger().log( std::to_string(source) + ": " + message + "\n" );
}


void callback_errorPrint(int error, const char* description)
{
    std::cout << description << " (error code " + std::to_string( error ) + ")\n";
}


void frameBufferResizeCallback( GLFWwindow* window,
                                int width,
                                int height  )
{
    std::cout << "Frame buffer resize!\n";
    glViewport( 0, 0, width, height );
}


} // namespace cie::gl