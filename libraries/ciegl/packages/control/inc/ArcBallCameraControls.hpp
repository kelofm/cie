#ifndef CIE_GL_ARCBALL_CAMERA_CONTROLS_HPP
#define CIE_GL_ARCBALL_CAMERA_CONTROLS_HPP

// --- Internal Includes ---
#include "packages/control/inc/MappedCameraControls.hpp"


namespace cie::gl {


/**
 * Configurable key controls:
 *  rotateUp ( "w", "up" )
 *  rotateDown ( "s", "down" )
 *  rotateLeft ( "a", "left" )
 *  rotateRight ( "d", "right" )
 *  moveForward ( "r" )
 *  moveBackward ( "f" )
 *  increaseZoomScale ( "kp_add" )
 *  decreaseZoomScale ( "kp_subtract" )
 */
class ArcBallCameraControls final : public MappedCameraControls
{
public:
    ArcBallCameraControls( bool verbose = false );

    const AbsCamera::vector_type& getCenter() const;
    void setCenter( const AbsCamera::vector_type& r_center );

    void setRotationScale( double scale );
    void setZoomScale( double scale );

private:
    void moveForward();
    void moveBackward();

    void rotate();
    void rotateLeft();
    void rotateRight();
    void rotateUp();
    void rotateDown();

    void zoomIn();
    void zoomOut();

    void increaseZoomScale();
    void decreaseZoomScale();

private:
    void onMouseButtonPress( KeyEnum button,
                             KeyEnum modifiers ) override;

    void onMouseButtonRelease( KeyEnum button,
                               KeyEnum modifiers ) override;

    void onCursorMovement( double x,
                           double y ) override;

    void onVerticalScroll( double offset ) override;

    ArcBallCameraControls::configuration_map_type makeConfigurationMap( ArcBallCameraControls::configuration_contents& r_configContents ) override;

private:
    AbsCamera::vector_type _center;

    double _rotationScale;
    double _zoomScale;
    double _zoomScaleIncrement;

    double _dx;
    double _dy;

    bool _isLeftMouseButtonPressed;
    bool _isRightMouseButtonPressed;

    bool _verbose;
};


} // namespace cie::gl


#endif