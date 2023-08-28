#ifndef CIE_GL_DEBUG_CAMERA_CONTROLS_HPP
#define CIE_GL_DEBUG_CAMERA_CONTROLS_HPP

// --- Internal Includes ---
#include "packages/control/inc/AbsCameraControls.hpp"


namespace cie::gl {


class DebugCameraControls final : public AbsCameraControls
{
public:
    DebugCameraControls();

    void onMouseButtonPress( KeyEnum button,
                             KeyEnum modifiers ) override;

    void onMouseButtonRelease( KeyEnum button,
                               KeyEnum modifiers ) override;

    void onCursorMovement( double x,
                           double y ) override;

    void onCursorEnter() override;

    void onCursorLeave() override;

    void onHorizontalScroll( double offset ) override;

    void onVerticalScroll( double offset ) override;

    void onKeyboardPress( KeyEnum key,
                          KeyEnum modifiers ) override;

    void onKeyboardRelease( KeyEnum key,
                            KeyEnum modifiers ) override;
};


} // namespace cie::gl


#endif