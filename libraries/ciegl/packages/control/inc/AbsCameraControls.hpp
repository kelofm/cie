#ifndef CIE_GL_ABS_CAMERA_CONTROLS_HPP
#define CIE_GL_ABS_CAMERA_CONTROLS_HPP

// --- Utility Includes ---
#include "packages/observer/inc/Observer.hpp"

// --- Internal Includes ---
#include "packages/context/inc/AbsWindow.hpp"
#include "packages/camera/inc/AbsCamera.hpp"
#include "packages/control/inc/callback_types.hpp"

// --- STL Includes ---
#include <functional>
#include <memory>
#include <set>


namespace cie::gl {


class AbsCameraControls : public utils::observer::Observer
{
public:
    AbsCameraControls();

    void mouseButtonCallback(KeyEnum button,
                             KeyEnum action,
                             KeyEnum modifiers);

    void cursorPositionCallback(double x,
                                double y);

    void cursorEnterCallback(KeyEnum entered);

    void scrollCallback(double xOffset,
                        double yOffset);

    void keyboardCallback(KeyEnum key,
                          KeyEnum scancode,
                          KeyEnum action,
                          KeyEnum modifiers);

    /// Bind window and camera to this instance
    void bind(AbsWindow::WeakPointer wp_window,
              CameraPtr p_camera);

protected:
    virtual void onMouseButtonPress([[maybe_unused]] KeyEnum button,
                                    [[maybe_unused]] KeyEnum modifiers)
    {}

    virtual void onMouseButtonHold([[maybe_unused]] KeyEnum button,
                                   [[maybe_unused]] KeyEnum modifiers)
    {}

    virtual void onMouseButtonRelease([[maybe_unused]] KeyEnum button,
                                      [[maybe_unused]] KeyEnum modifiers)
    {}

    virtual void onCursorMovement([[maybe_unused]] double x,
                                  [[maybe_unused]] double y)
    {}

    virtual void onCursorEnter()
    {}

    virtual void onCursorLeave()
    {}

    virtual void onHorizontalScroll([[maybe_unused]] double offset)
    {}

    virtual void onVerticalScroll([[maybe_unused]] double offset)
    {}

    virtual void onKeyboardPress(KeyEnum key,
                                  KeyEnum modifiers);

    virtual void onKeyboardHold([[maybe_unused]] KeyEnum key,
                                [[maybe_unused]] KeyEnum modifiers)
    {}

    virtual void onKeyboardRelease([[maybe_unused]] KeyEnum key,
                                   [[maybe_unused]] KeyEnum modifiers)
    {}

    virtual void onSubjectChange() override
    {}

protected:
    AbsWindow::SharedPointer getWindow();

protected:
    AbsWindow::WeakPointer _p_window;

    CameraPtr _p_camera;

    double _x;

    double _y;

    KeyEnum _escapeKey;

    std::set<KeyEnum> _activeKeys;
};


using CameraControlsPtr = std::shared_ptr<AbsCameraControls>;


} // namespace cie::gl


#endif
