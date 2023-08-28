// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/context/inc/GLFWMonitor.hpp"


namespace cie::gl {


GLFWMonitor::GLFWMonitor( Size id,
                          Size width,
                          Size height,
                          bool isPrimary,
                          GLFWmonitor* p_monitor ) :
    AbsMonitor( id,
                width,
                height,
                isPrimary ),
    _p_monitor( p_monitor )
{
}


GLFWMonitor::~GLFWMonitor()
{
}


MonitorPtr makePrimaryGLFWMonitor()
{
    CIE_BEGIN_EXCEPTION_TRACING

    GLFWmonitor* p_monitor = glfwGetPrimaryMonitor();
    if ( !p_monitor )
        CIE_THROW( Exception, "Failed to retrieve primary monitor" )

    const auto* p_videoMode = glfwGetVideoMode( p_monitor );

    return MonitorPtr(
        new GLFWMonitor( 0,
                         Size( p_videoMode->width ),
                         Size( p_videoMode->height ),
                         true,
                         p_monitor )
    );

    CIE_END_EXCEPTION_TRACING
}


std::shared_ptr<std::vector<MonitorPtr>> makeAllGLFWMonitors()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_monitors = std::make_shared<std::vector<MonitorPtr>>();
    p_monitors->push_back( makePrimaryGLFWMonitor() );

    int numberOfMonitors            = 0;
    GLFWmonitor** pp_monitor        = glfwGetMonitors( &numberOfMonitors );
    GLFWmonitor* p_primaryMonitor   = glfwGetPrimaryMonitor();

    p_monitors->reserve( numberOfMonitors );

    for ( Size i=0; i<Size(numberOfMonitors); ++i )
    {
        GLFWmonitor* p_currentMonitor = *(pp_monitor + i);

        if ( p_currentMonitor == p_primaryMonitor )
            continue;

        const auto* p_videoMode = glfwGetVideoMode( p_currentMonitor );
        p_monitors->push_back( MonitorPtr(
            new GLFWMonitor( i+1,
                             Size( p_videoMode->width ),
                             Size( p_videoMode->height ),
                             false,
                             p_currentMonitor )
        ));
    }

    return p_monitors;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl