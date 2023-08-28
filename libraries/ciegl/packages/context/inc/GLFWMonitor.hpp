#ifndef CIE_GL_GLFW_MONITOR_HPP
#define CIE_GL_GLFW_MONITOR_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- Intenral Includes ---
#include "packages/context/inc/AbsMonitor.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>


namespace cie::gl {


class GLFWMonitor final : public AbsMonitor
{
public:
    ~GLFWMonitor() override;

    friend MonitorPtr makePrimaryGLFWMonitor();
    friend std::shared_ptr<std::vector<MonitorPtr>> makeAllGLFWMonitors();

private:
    GLFWMonitor( Size id,
                 Size width,
                 Size height,
                 bool isPrimary,
                 GLFWmonitor* p_monitor );

private:
    GLFWmonitor* _p_monitor;
};


} // namespace cie::gl


#endif