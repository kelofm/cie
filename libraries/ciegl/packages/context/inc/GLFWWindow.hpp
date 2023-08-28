#ifndef CIE_GL_GLFW_WINDOW_HPP
#define CIE_GL_GLFW_WINDOW_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- Internal Includes ---
#include "packages/context/inc/AbsWindow.hpp"
#include "packages/control/inc/callback_types.hpp"


namespace cie::gl {


class GLFWWindow final : public AbsWindow
{
public:
    friend class GLFWContext;

public:
    ~GLFWWindow();

    virtual void activate() override;

    const GLFWwindow* get() const;

    GLFWwindow* get();

protected:
    virtual void setSizeImpl(WindowSize size) override;

    virtual bool updateImpl() override;

    virtual void setWindowName(std::string::const_iterator begin,
                               std::string::const_iterator end) override;

private:
    GLFWWindow(GLFWwindow* p_glfwWindow,
               Size id,
               const std::string& r_name,
               Size width,
               Size height,
               utils::Logger& r_logger = utils::LoggerSingleton::get());

    friend void GLFWMouseButtonCallback(GLFWwindow* p_window, KeyEnum button, KeyEnum action, KeyEnum modifier);

    friend void GLFWCursorPositionCallback(GLFWwindow* p_window, Double x, Double y);

    friend void GLFWCursorEnterCallback(GLFWwindow* p_window, KeyEnum hasEntered);

    friend void GLFWScrollCallback(GLFWwindow* p_window, Double xOffset, Double yOffset);

    friend void GLFWKeyboardCallback(GLFWwindow* p_window, KeyEnum key, KeyEnum scancode, KeyEnum action, KeyEnum modifier);

    friend void GLFWWindowResizeCallback(GLFWwindow* p_window, int width, int height);

private:
    GLFWwindow* _p_window;
};


namespace detail {
/** Helper function for getting a GLFW window
 * from a AbsWindow::SharedPointer interface
*/
GLFWwindow* getGLFWwindow(AbsWindow::SharedPointer p_window);
} // namespace detail


} // namespace cie::gl

#include "packages/context/impl/GLFWWindow_impl.hpp"

#endif