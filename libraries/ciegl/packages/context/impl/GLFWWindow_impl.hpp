#ifndef CIE_GL_GLFW_WINDOW_IMPL_HPP
#define CIE_GL_GLFW_WINDOW_IMPL_HPP

// --- GL Includes ---
#include "packages/context/inc/GLFWWindow.hpp"


namespace cie::gl {


inline void GLFWMouseButtonCallback(GLFWwindow* p_window, KeyEnum button, KeyEnum action, KeyEnum modifiers)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerMouseButtonCallback(button, action, modifiers);
}


inline void GLFWCursorPositionCallback(GLFWwindow* p_window, Double x, Double y)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerCursorPositionCallback(x, y);
}


inline void GLFWCursorEnterCallback(GLFWwindow* p_window, KeyEnum hasEntered)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerCursorEnterCallback(hasEntered);
}


inline void GLFWScrollCallback(GLFWwindow* p_window, Double xOffset, Double yOffset)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerScrollCallback(xOffset, yOffset);
};


inline void GLFWKeyboardCallback(GLFWwindow* p_window, KeyEnum key, KeyEnum scancode, KeyEnum action, KeyEnum modifier)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerKeyboardCallback(key, scancode, action, modifier);
};


inline void GLFWWindowResizeCallback(GLFWwindow* p_window, int width, int height)
{
    auto p_this = static_cast<GLFWWindow*>(glfwGetWindowUserPointer(p_window));
    CIE_CHECK_POINTER(p_this)
    p_this->triggerWindowResizeCallback(width, height);
};


} // namespace cie::gl


#endif