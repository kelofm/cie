// --- External Includes ---
#include "glad/glad.h"
#include "GLFW/glfw3.h"

// --- Internal Includes ---
#include "packages/control/inc/DebugCameraControls.hpp"


namespace cie::gl {


DebugCameraControls::DebugCameraControls() :
    AbsCameraControls()
{
}


void DebugCameraControls::onMouseButtonPress( KeyEnum button,
                                              KeyEnum modifiers )
{
    this->_p_camera->log(
        "Mouse press '" + std::to_string(button) + "' with modifiers: '" + std::to_string(modifiers) + "'"
    );
}


void DebugCameraControls::onMouseButtonRelease( KeyEnum button,
                                                KeyEnum modifiers )
{
    this->_p_camera->log(
        "Mouse release '" + std::to_string(button) + "' with modifiers: '" + std::to_string(modifiers) + "'"
    );
}


void DebugCameraControls::onCursorMovement( double x,
                                            double y )
{
    this->_p_camera->log(
        "Mouse position :" + std::to_string(x) + ", " + std::to_string(y)
    );
}


void DebugCameraControls::onCursorEnter()
{
    this->_p_camera->log( "Cursor enters window" );
}


void DebugCameraControls::onCursorLeave()
{
    this->_p_camera->log( "Cursor leaves window" );
}


void DebugCameraControls::onHorizontalScroll( double offset )
{
    this->_p_camera->log( "Horizontal scroll: " + std::to_string(offset) );
}


void DebugCameraControls::onVerticalScroll( double offset )
{
    this->_p_camera->log( "Vertical scroll: " + std::to_string(offset) );
}


void DebugCameraControls::onKeyboardPress( KeyEnum key,
                                           KeyEnum modifiers )
{
    this->_p_camera->log(
        "Press key '" + std::to_string(key) + "' with modifiers: '" + std::to_string(modifiers) + "'"
    );

    if ( key == GLFW_KEY_ESCAPE )
    {
        auto p_window = this->getWindow();
        p_window->endLoop();
    }
}


void DebugCameraControls::onKeyboardRelease( KeyEnum key,
                                             KeyEnum modifiers )
{
    this->_p_camera->log(
        "Release key '" + std::to_string(key) + "' with modifiers: '" + std::to_string(modifiers) + "'"
    );
}


} // namespace cie::gl