#ifndef CIE_GL_FLY_CAMERA_CONTROLS_HPP
#define CIE_GL_FLY_CAMERA_CONTROLS_HPP

// --- Internal Includes ---
#include "packages/control/inc/MappedCameraControls.hpp"


namespace cie::gl {


/**
 * Configurable key controls:
 *  moveForward ( "w", "up" )
 *  moveBackward ( "s", "down" )
 *  strafeLeft ( "a", "left" )
 *  strafeRight ( "d", "right" )
 *  moveUp ( "r", "space" )
 *  moveDown ( "f", "left_control" )
 *  rollCounterClockwise( "q" )
 *  rollClockwise( "r" )
 */
class FlyCameraControls final : public MappedCameraControls
{
public:
    FlyCameraControls( bool verbose = false );

    void setMovementScale( double scale );
    void setRotationScale( double scale );
    void setZoomScale( double scale );

private:
    void moveForward();
    void moveBackward();
    void strafeLeft();
    void strafeRight();
    void moveUp();
    void moveDown();

    void rotate();
    void rollCounterClockwise();
    void rollClockwise();

    void zoomIn();
    void zoomOut();

private:
    void onMouseButtonPress( KeyEnum button,
                             KeyEnum modifiers ) override;

    void onMouseButtonRelease( KeyEnum button,
                               KeyEnum modifiers ) override;

    void onCursorMovement( double x,
                           double y ) override;

    void onHorizontalScroll( double offset ) override;

    void onVerticalScroll( double offset ) override;

    FlyCameraControls::configuration_map_type makeConfigurationMap( FlyCameraControls::configuration_contents& r_configContents ) override;

private:
    double _movementScale;
    double _rotationScale;
    double _zoomScale;

    double _x;
    double _y;
    double _dx;
    double _dy;

    bool _isLeftMouseButtonPressed;
    bool _isRightMouseButtonPressed;

    bool _verbose;
};


} // namespace cie::gl


#endif