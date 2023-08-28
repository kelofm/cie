// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/context/inc/defaultCallback.hpp"

// --- STL Includes ---
#include <string>
#include <iostream>

namespace cie::gl {


/* --- GLFWContext --- */

GLFWContext::GLFWContext(Size versionMajor,
                         Size versionMinor,
                         Size MSAASamples,
                         utils::LoggerPtr p_logger)   :
    AbsContext(versionMajor,
                versionMinor,
                MSAASamples,
                p_logger),
    _isGLADInitialized(false),
    _isCurrent(false),
    _mutex()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto scopedBlock = this->logger().newBlock("GLFWContext init");

    Lock lock(_mutex);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, versionMajor);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, versionMinor);
    glfwWindowHint(GLFW_SAMPLES, MSAASamples);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifndef NDEBUG
        glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
    #endif

    glfwSetErrorCallback(callback_errorPrint);

    if (!glfwInit())
        this->logger().error<Exception>("Failed to initialize GLFW");
    else
        this->logger().log("GLFW initialization successful");

    this->logger().log("Open context");

    CIE_END_EXCEPTION_TRACING
}


GLFWContext::~GLFWContext()
{
    glfwTerminate();
}



AbsWindow::SharedPointer GLFWContext::newWindowImpl(Size width,
                                                    Size height,
                                                    const std::string& r_name)
{
    CIE_BEGIN_EXCEPTION_TRACING

    GLFWwindow* p_glfwWindow;
    {
        Lock lock(_mutex);
        p_glfwWindow = glfwCreateWindow(width, height, r_name.c_str(), nullptr, nullptr);
        CIE_CHECK(p_glfwWindow, "Failed to create window " << r_name)
    }
    this->makeCurrent(p_glfwWindow);

    auto p_window = AbsWindow::SharedPointer(
        new GLFWWindow(p_glfwWindow,
                       0,
                       r_name,
                       width,
                       height,
                       this->logger())
    );
    return p_window;

    CIE_END_EXCEPTION_TRACING
}



void GLFWContext::focusWindow(AbsWindow::SharedPointer p_window)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (p_window != nullptr) [[likely]]
        this->logger().log("Focus on window '" + p_window->name() + "'");
    else
        this->logger().error<Exception>("Focus on non-existent window!");

    CIE_END_EXCEPTION_TRACING
}


void GLFWContext::makeCurrent(GLFWwindow* p_window)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto block = this->logger().newBlock("Make context current");

    CIE_CHECK_POINTER(p_window)

    if (!_isCurrent)
    {
        {
            Lock lock(_mutex);
            glfwMakeContextCurrent(p_window);
            _isCurrent = true;
        }
        this->logger().log("Context is now current");
    }
    else
        this->logger().warn("Making current context current again");

    this->ensureGLADInitialized();

    {
        Lock lock(_mutex);
        glDebugMessageCallback(messageCallback, this);
    }

    CIE_END_EXCEPTION_TRACING
}



void GLFWContext::closeWindowImpl(AbsWindow::SharedPointer p_window)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->logger().log("Force close window " + p_window->name());
    Lock lock(_mutex);
    glfwSetWindowShouldClose(detail::getGLFWwindow(p_window), 1);
    CIE_END_EXCEPTION_TRACING
}


void GLFWContext::ensureGLADInitialized()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto block = this->logger().newBlock("GLAD initialization");

    Lock lock(_mutex);
    if (_isGLADInitialized)
    {
        this->logger().log("GLAD is already initialized");
        return;
    }

    if ((!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)))
        this->logger().error<Exception>("Failed to initialize GLAD");

    _isGLADInitialized = true;

    CIE_END_EXCEPTION_TRACING
}




/* --- GLFWContextSingleton --- */

GLFWContext::SharedPointer GLFWContextSingleton::_p_context = nullptr;


std::mutex GLFWContextSingleton::_mutex;


GLFWContext::WeakPointer GLFWContextSingleton::get(Size versionMajor,
                                                   Size versionMinor,
                                                   Size MSAASamples,
                                                   utils::LoggerPtr p_logger)
{
    CIE_BEGIN_EXCEPTION_TRACING

    Lock lock(_mutex);

    if (GLFWContextSingleton::_p_context) // context is already initialized
    {
        // Check context parameters
        if (
            GLFWContextSingleton::_p_context->version().first != versionMajor
            ||
            GLFWContextSingleton::_p_context->version().second != versionMinor
            ||
            GLFWContextSingleton::_p_context->MSAASamples() != MSAASamples
        )
            GLFWContextSingleton::_p_context->logger().error("GLFWContext parameter mismatch!");

        // TODO: Add stream if necessary
    }
    else // initialize the context
    {
        GLFWContextSingleton::_p_context = GLFWContext::SharedPointer(new GLFWContext(
            versionMajor,
            versionMinor,
            MSAASamples,
            p_logger
        ));
    }

    return GLFWContextSingleton::_p_context;

    CIE_END_EXCEPTION_TRACING
}


GLFWContext::WeakPointer GLFWContextSingleton::get()
{
    CIE_BEGIN_EXCEPTION_TRACING

    GLFWContext::SharedPointer p_context;
    //Lock lock(_mutex);

    if (GLFWContextSingleton::_p_context)
        p_context = GLFWContextSingleton::_p_context;
    else
        p_context = GLFWContextSingleton::get(4, 5, 0, utils::LoggerSingleton::getPtr()).lock();

    return p_context;

    CIE_END_EXCEPTION_TRACING
}


GLFWContext::WeakPointer GLFWContextSingleton::get(utils::LoggerPtr p_logger)
{
    CIE_BEGIN_EXCEPTION_TRACING

    GLFWContext::SharedPointer p_context;
    //Lock lock(_mutex);

    if (!GLFWContextSingleton::_p_context)
        p_context = GLFWContextSingleton::get(4, 5, 0, p_logger).lock();
    else
    {
        if (&GLFWContextSingleton::_p_context->logger() != &*p_logger)
            CIE_THROW(Exception, "Cannot assign new logger to GLFWContextSingleton")

        p_context = GLFWContextSingleton::_p_context;
    }

    return p_context;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl