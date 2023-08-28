// --- Internal Includes ---
#include "packages/context/inc/windowutilities.hpp"

// --- STL Includes ---
#include <utility>

namespace cie::gl {


std::pair<size_t, size_t> getMonitorResolution( GLFWmonitor* monitor )
{
    const GLFWvidmode* mode = glfwGetVideoMode( monitor );
    return std::make_pair<size_t, size_t>(size_t(mode->width), size_t(mode->height));
}


std::pair<size_t, size_t> getPrimaryMonitorResolution()
{
    const GLFWvidmode* mode = glfwGetVideoMode( glfwGetPrimaryMonitor() );
    return std::make_pair<size_t, size_t>( size_t(mode->width), size_t(mode->height) );
}


} // namespace cie::gl