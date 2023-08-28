// --- External Includes ---
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/control/inc/AbsCameraControls.hpp"

// --- STL Includes ---
#include <functional>


namespace cie::gl {


AbsCameraControls::AbsCameraControls() :
    utils::observer::Observer(),
    _p_window(),
    _p_camera(nullptr),
    _x(0.0),
    _y(0.0),
    _escapeKey(GLFW_KEY_ESCAPE),
    _activeKeys()
{
}


void AbsCameraControls::mouseButtonCallback(KeyEnum button,
                                             KeyEnum action,
                                             KeyEnum modifiers)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(this->_p_camera)

    if (action == GLFW_PRESS || action == GLFW_REPEAT)
        this->onMouseButtonPress(button, modifiers);
    else if (action == GLFW_RELEASE)
        this->onMouseButtonRelease(button, modifiers);
    else
        _p_camera->log("Unrecognized mouse button action: " + std::to_string(action),
                        LOG_TYPE_ERROR);

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::cursorPositionCallback(double x,
                                                double y)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(this->_p_camera)

    this->onCursorMovement(x, y);

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::cursorEnterCallback(KeyEnum entered)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(this->_p_camera)

    if (entered == GLFW_TRUE)
        this->onCursorEnter();
    else if (entered == GLFW_FALSE)
        this->onCursorLeave();
    else
        _p_camera->log("Unrecognised cursor enter action: " + std::to_string(entered),
                        LOG_TYPE_ERROR);

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::scrollCallback(double xOffset,
                                        double yOffset)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(this->_p_camera)

    if (xOffset != 0.0)
        this->onHorizontalScroll(xOffset);
    if (yOffset != 0.0)
        this->onVerticalScroll(yOffset);

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::keyboardCallback(KeyEnum key,
                                         KeyEnum scancode,
                                         KeyEnum action,
                                         KeyEnum modifiers)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(this->_p_camera)

    if (key == this->_escapeKey)
    {
        this->getWindow()->endLoop();
        return;
    }

    if (action==GLFW_PRESS || action==GLFW_REPEAT)
    {
        this->_activeKeys.insert(key);
        this->onKeyboardPress(key, modifiers);
    }
    else if (action == GLFW_RELEASE)
    {
        this->_activeKeys.erase(key);
        this->onKeyboardRelease(key, modifiers);
    }
    else
        this->_p_camera->log("Unrecognised keyboard action: " + std::to_string(action),
                              LOG_TYPE_ERROR);

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::onKeyboardPress(KeyEnum key,
                                         KeyEnum modifiers)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (key == GLFW_KEY_ESCAPE)
    {
        auto p_window = this->getWindow();
        p_window->endLoop();
    }

    CIE_END_EXCEPTION_TRACING
}


void AbsCameraControls::bind(AbsWindow::WeakPointer wp_window,
                              CameraPtr p_camera)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (auto p_window = wp_window.lock())
    {
        CIE_CHECK_POINTER(p_camera)

        this->_p_window = p_window;
        this->_p_camera = p_camera;

        p_window->log("Bind camera controls");
        this->setSubject(p_window);
        p_window->attachObserver(this->shared_from_this());
        p_window->setMouseButtonCallback(std::bind(&AbsCameraControls::mouseButtonCallback, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
        p_window->setCursorPositionCallback(std::bind(&AbsCameraControls::cursorPositionCallback, this, std::placeholders::_1, std::placeholders::_2));
        p_window->setCursorEnterCallback(std::bind(&AbsCameraControls::cursorEnterCallback, this, std::placeholders::_1));
        p_window->setScrollCallback(std::bind(&AbsCameraControls::scrollCallback, this, std::placeholders::_1, std::placeholders::_2));
        p_window->setKeyboardCallback([this](KeyEnum key, KeyEnum scancode, KeyEnum action, KeyEnum modifier) {this->keyboardCallback(key, scancode, action, modifier);});
    }
    else
        CIE_THROW(Exception, "Invalid window pointer")

    CIE_END_EXCEPTION_TRACING
}


AbsWindow::SharedPointer AbsCameraControls::getWindow()
{
    if (auto p_window = this->_p_window.lock())
        return p_window;
    else
        CIE_THROW(Exception, "Pointer to window expired!")

    return nullptr;
}


} // namespace cie::gl