// --- Utility Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- GL Includes ---
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/utility/inc/GLError.hpp"
#include "packages/context/inc/GLFWContext.hpp"


namespace cie::gl {


GLFWWindow::GLFWWindow(GLFWwindow* p_glfwWindow,
                       Size id,
                       const std::string& r_name,
                       Size width,
                       Size height,
                       utils::Logger& r_logger) :
    AbsWindow(id,
              r_name,
              width,
              height,
              r_logger),
              _p_window(p_glfwWindow)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto block = this->newBlock("Constructing window '" + r_name + "'");

    // Set window pointer
    glfwSetWindowUserPointer(_p_window, this);
    checkGLErrors("Failed to set window user pointer for '" + r_name + "'");

    // Set callback rerouters
    glfwSetMouseButtonCallback(_p_window, GLFWMouseButtonCallback);
    glfwSetCursorPosCallback(_p_window, GLFWCursorPositionCallback);
    glfwSetCursorEnterCallback(_p_window, GLFWCursorEnterCallback);
    glfwSetScrollCallback(_p_window, GLFWScrollCallback);
    glfwSetKeyCallback(_p_window, GLFWKeyboardCallback);
    glfwSetWindowSizeCallback(_p_window, GLFWWindowResizeCallback);
    checkGLErrors("Failed to set callbacks for window '" + r_name + "'");

    this->setSize(this->getSize());
    this->setName(r_name);

    // Check whether the new window has the requested parameters
    int checkWidth, checkHeight;
    glfwGetFramebufferSize(_p_window, &checkWidth, &checkHeight);
    if (Size(checkWidth)!=width || Size(checkHeight)!=height)
        CIE_THROW(Exception, "Created window is not of the requested size!");

    checkGLErrors("Failed to construct window '" + r_name + "'");

    CIE_END_EXCEPTION_TRACING
}


GLFWWindow::~GLFWWindow()
{
    glfwDestroyWindow(this->_p_window);
}


void GLFWWindow::activate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (auto p_context = GLFWContextSingleton::get().lock())
        p_context->makeCurrent(_p_window);
    else
        CIE_THROW(NullPtrException, "")

    CIE_END_EXCEPTION_TRACING
}


void GLFWWindow::setSizeImpl(WindowSize size)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(_p_window)
    glfwSetWindowSize(_p_window, size.first, size.second);

    CIE_END_EXCEPTION_TRACING
}


bool GLFWWindow::updateImpl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check for external errors
    checkGLErrors(*this,
                   "Error before updating!");

    glfwPollEvents();

    // Check for polling errors
    checkGLErrors(*this,
                   "Error polling events!");

    glfwSwapBuffers(this->_p_window);

    // Check for buffer errors
    checkGLErrors(*this,
                   "Error swapping buffers!");

    return true;

    CIE_END_EXCEPTION_TRACING
}


void GLFWWindow::setWindowName(std::string::const_iterator begin,
                               std::string::const_iterator end)
{
    glfwSetWindowTitle(_p_window, std::string(begin, end).c_str());
}


const GLFWwindow* GLFWWindow::get() const
{
    return _p_window;
}


GLFWwindow* GLFWWindow::get()
{
    return _p_window;
}


namespace detail {
GLFWwindow* getGLFWwindow(AbsWindow::SharedPointer p_window)
{
    return dynamic_pointer_cast<GLFWWindow>(p_window)->get();
}
} // namespace detail


} // namespace cie::gl