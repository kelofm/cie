#ifndef CIE_GL_PAN_ZOOM_CAMERA_CONTROLS_HPP
#define CIE_GL_PAN_ZOOM_CAMERA_CONTROLS_HPP

// --- Internal Includes ---
#include "packages/control/inc/MappedCameraControls.hpp"


namespace cie::gl {


/**
 * Configurable key controls:
 *  moveUp ( "w", "up" )
 *  moveDown ( "s", "down" )
 *  moveLeft ( "a", "left" )
 *  moveRight ( "d", "right" )
 *  zoomIn ( "kp_add" )
 *  zoomOut ( "kp_subtract" )
 */
class PanZoomCameraControls final : public MappedCameraControls
{
public:
    PanZoomCameraControls( bool verbose = false );

    void setMovementScale( double scale );
    void setZoomScale( double scale );

private:
    void moveHorizontally( double amount );
    void moveVertically( double amount );

    void moveUp();
    void moveDown();
    void moveLeft();
    void moveRight();

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

    PanZoomCameraControls::configuration_map_type makeConfigurationMap( PanZoomCameraControls::configuration_contents& r_configContents ) override;

private:
    double _movementScale;
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