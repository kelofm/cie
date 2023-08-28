#ifndef CIEGL_WINDOW_UTILITIES_HPP
#define CIEGL_WINDOW_UTILITIES_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- STL Includes ---
#include <utility>
#include <cstddef>

namespace cie::gl {


std::pair<size_t, size_t> getMonitorResolution( GLFWmonitor* monitor=glfwGetPrimaryMonitor() );
std::pair<size_t, size_t> getPrimaryMonitorResolution();


} // namespace cie::gl


#endif