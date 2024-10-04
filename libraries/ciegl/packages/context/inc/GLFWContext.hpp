#ifndef CIEGL_GL_GLFW_CONTEXT_HPP
#define CIEGL_GL_GLFW_CONTEXT_HPP

// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- Utility Includes ---
#include "packages/logging/inc/Logger.hpp"
#include "packages/observer/inc/Observer.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- Internal Includes ---
#include "packages/context/inc/AbsContext.hpp"
#include "packages/context/inc/windowutilities.hpp"
#include "packages/context/inc/defaultCallback.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <functional>
#include <memory>
#include <string>
#include <filesystem>
#include <mutex>


namespace cie::gl {


/// @brief Class representing a GLFW context.
///
/// @note Only to be used with @ref GLFWContextSingleton.
/// @ingroup ciegl
class GLFWContext final : public AbsContext
{
private:
    friend class GLFWContextSingleton;

    using Mutex = std::mutex;

    using Lock = std::scoped_lock<Mutex>;

    CIE_DEFINE_CLASS_POINTERS(GLFWContext)

public:
    ~GLFWContext();

    /// @brief Activate the specified window.
    ///
    /// @warning The window must exist, but this is not checked.
    /// @todo Implement and test.
    void focusWindow(AbsWindow::SharedPointer p_window) override;

    /** @brief Make this context the current one.
     *
     *  @warning Must be called @b after @p glfwCreateWindow but @b before any other GL calls,
     *           which includes the construction of a @p GLFWWindow. This is why the constructor
     *           of @p GLFWWindow is private.
     *  @todo Eliminate repeated calls, and check for bugs.
     */
    void makeCurrent(GLFWwindow* p_window);

private:
    /** @brief Initialize the GLFW context.
     *
     *  @param versionMajor: OpenGL standard version major specifier (must match the version of the existing context).
     *  @param versionMinor: OpenGL standard version minor specifier (must match the version of the existing context).
     *  @param MSAASamples: number of MSAASamples (must match the version of the existing context).
     *  @param p_logger: Existing logger to use. Default gets the static logger from @utils::LoggerSingleton.
     *  @note This constructor is private and should not be called directly from anywhere
     *        except from @ref GLFWContextSingleton. See the overloads of @ref GLFWContextSingleton::get.
     */
    GLFWContext(Size versionMajor = 4,
                Size versionMinor = 5,
                Size MSAASamples = 0,
                utils::LoggerPtr p_logger = utils::LoggerSingleton::getPtr());

    /// Check whether GLAD was initialized and initialize it if it wasn't
    void ensureGLADInitialized();

    AbsWindow::SharedPointer newWindowImpl(size_t width = 800,
                                           size_t height = 600,
                                           const std::string& r_name = "CiE") override;

    void closeWindowImpl(AbsWindow::SharedPointer p_window) override;

private:
    bool _isGLADInitialized;

    bool _isCurrent;

    Mutex _mutex;
};



/** @brief GLFWContext manager ensuring a single context at all times.
 *
 *  @bug Something isn't right between created windows.
 *  @ingroup ciegl
 */
class GLFWContextSingleton
{
private:
    using Mutex = GLFWContext::Mutex;

    using Lock = GLFWContext::Lock;

public:
    /** @brief Get the existing context or construct a new one.
     *
     *  @details Check whether the GLFWContext was initialized with the specified parameters and create one if it was not.
     *  @param versionMajor: OpenGL standard version major specifier (must match the version of the existing context).
     *  @param versionMinor: OpenGL standard version minor specifier (must match the version of the existing context).
     *  @param MSAASamples: number of MSAASamples (must match the version of the existing context).
     *  @param r_logFilePath: Path to the log file.
     *  @param useConsole: Print log output to std::cout.
     *
     *  @return A pointer to the static GLFW context.
     *  @note Throws an error if the requested parameters do not match the
     *        parameters of the existing context.
     */
    static GLFWContext::WeakPointer get(Size versionMajor,
                                        Size versionMinor,
                                        Size MSAASamples,
                                        utils::LoggerPtr p_logger = utils::LoggerSingleton::getPtr());

    /// @brief Get the existing context or create a default one.
    static GLFWContext::WeakPointer get();

    /** @brief Initialize the GLFW context or return the existing one, and add a log file to its list of streams.
     *
     *  @param r_logFilePath: path to the log file.
     *  @return a pointer to the GLFW context.
     *  @todo untested.
     */
    static GLFWContext::WeakPointer get(utils::LoggerPtr p_logger);

private:
    static GLFWContext::SharedPointer _p_context;

    static std::mutex _mutex;
};



} // namespace cie::gl

#endif