// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/context/inc/GLFWMonitor.hpp"
#include "cmake_variables.hpp"


namespace cie::gl {


CIE_TEST_CASE( "GLFWContext", "[context]" )
{
    CIE_TEST_CASE_INIT( "GLFWContext" )

    // Window parameters
    std::pair<Size,Size> windowSize   { 800, 600 };
    const std::string windowName    = "GLFW testwindow";

    {
        CIE_TEST_CASE_INIT( "Context creation" )

        ContextPtr p_context;

        CIE_TEST_REQUIRE_NOTHROW(
            p_context = gl::GLFWContextSingleton::get().lock()
        );
    }

    {
        CIE_TEST_CASE_INIT( "Window creation" )

        auto p_context = gl::GLFWContextSingleton::get();

        AbsWindow::SharedPointer p_window;
        CIE_TEST_REQUIRE_NOTHROW(
            p_window = p_context.lock()->newWindow(
                windowSize.first,
                windowSize.second,
                windowName
            )
        );

        CIE_TEST_CHECK_NOTHROW( p_context.lock()->focusWindow(p_window) );
    }

}


} // namespace cie::gl