#ifndef CIEGL_CALLBACKS_KEYS_HPP
#define CIEGL_CALLBACKS_KEYS_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>


namespace cie::gl {


void GLAPIENTRY messageCallback( GLenum source,
                                 GLenum type,
                                 GLuint id,
                                 GLenum severity,
                                 GLsizei length,
                                 const GLchar* message,
                                 const void* logger );
                            
void callback_errorPrint(int error, const char* description);

void frameBufferResizeCallback( GLFWwindow* window,
                                int width,
                                int height  );


} // namespace cie::gl

#endif