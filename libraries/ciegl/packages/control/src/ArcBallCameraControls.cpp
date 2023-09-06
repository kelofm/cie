#define _USE_MATH_DEFINES

// --- External Includes ---
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/control/inc/ArcBallCameraControls.hpp"
#include "packages/utility/inc/glm_overloads.hpp"
#include "packages/control/inc/utilities.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <math.h>


namespace cie::gl {


ArcBallCameraControls::ArcBallCameraControls( bool verbose ) :
    MappedCameraControls(),
    _center( {0.0, 0.0, 0.0} ),
    _rotationScale( 0.1 * M_PI / 180.0 ),
    _zoomScale( 1.05 ),
    _zoomScaleIncrement( 0.05 ),
    _dx( 0.0 ),
    _dy( 0.0 ),
    _isLeftMouseButtonPressed( false ),
    _isRightMouseButtonPressed( false ),
    _verbose( verbose )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Set default controls
    std::filesystem::path confiFilePath = getDataPath() / "ciegl/configuration/controls/ArcBallCameraControls.json";
    std::ifstream configFile( confiFilePath );

    if ( !configFile.is_open() )
        CIE_THROW( Exception, "Failed to open control config file: " + confiFilePath.string() );

    this->configure( configFile );

    CIE_END_EXCEPTION_TRACING
}


const AbsCamera::vector_type& ArcBallCameraControls::getCenter() const
{
    return this->_center;
}


void ArcBallCameraControls::setCenter( const AbsCamera::vector_type& r_center )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate( r_center - this->_center );
    this->_center = r_center;

    if ( this->_verbose )
        *this->_p_camera  << "Set center to" << r_center;

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::setRotationScale( double scale )
{
    if ( scale < 0.0 )
        this->_p_camera->log( "Rotation scale must be non-negative (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );

    this->_rotationScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set rotation scale to " + std::to_string( scale );
}


void ArcBallCameraControls::setZoomScale( double scale )
{
    if ( scale <= 0.0 )
        this->_p_camera->log( "Zoom scale must be positive (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );
    this->_zoomScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set zoom scale to " + std::to_string( scale );
}


void ArcBallCameraControls::moveForward()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        this->_zoomScale * this->_p_camera->direction()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move forward to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::moveBackward()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        -this->_zoomScale * this->_p_camera->direction()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move backward to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::rotate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Rotate 'horizontally'
    this->_p_camera->rotate(
        -this->_rotationScale * this->_dx,
        { 0.0, 0.0, 1.0 },
        this->_center
    );

    // Rotate 'vertically'
    auto axis = glm::cross(
        { 0.0, 0.0, 1.0 },
        this->_p_camera->direction()
    );

    this->_p_camera->rotate(
        -this->_rotationScale * this->_dy,
        axis,
        this->_center
    );

    if ( this->_verbose )
        *this->_p_camera << "Rotate to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::rotateLeft()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->rotate(
        -this->_rotationScale,
        { 0.0, 0.0, 1.0 },
        this->_center
    );

    if ( this->_verbose )
        *this->_p_camera << "Rotate left to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::rotateRight()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->rotate(
        this->_rotationScale,
        { 0.0, 0.0, 1.0 },
        this->_center
    );

    if ( this->_verbose )
        *this->_p_camera << "Rotate right to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::rotateUp()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto axis = glm::cross(
        { 0.0, 0.0, 1.0 },
        this->_p_camera->direction()
    );

    this->_p_camera->rotate(
        this->_rotationScale,
        axis,
        this->_center
    );

    if ( this->_verbose )
        *this->_p_camera << "Rotate up to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::rotateDown()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto axis = glm::cross(
        { 0.0, 0.0, 1.0 },
        this->_p_camera->direction()
    );

    this->_p_camera->rotate(
        -this->_rotationScale,
        axis,
        this->_center
    );

    if ( this->_verbose )
        *this->_p_camera << "Rotate down to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::zoomIn()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom in";

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::zoomOut()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( 1.0 / this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom out";

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::increaseZoomScale()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_zoomScale += this->_zoomScaleIncrement;

    if ( this->_verbose )
        *this->_p_camera << "Increase zoom scale to " + std::to_string( this->_zoomScale );

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::decreaseZoomScale()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_zoomScale -= this->_zoomScaleIncrement;

    if ( this->_zoomScale <= 1.0 )
        this->_zoomScale = 1.0 + this->_zoomScaleIncrement;

    if ( this->_verbose )
        *this->_p_camera << "Decrease zoom scale to " + std::to_string( this->_zoomScale );

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::onMouseButtonPress(KeyEnum button,
                                               [[maybe_unused]] KeyEnum modifiers)
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = true;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = true;
}


void ArcBallCameraControls::onMouseButtonRelease(KeyEnum button,
                                                 [[maybe_unused]] KeyEnum modifiers)
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = false;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = false;
}


void ArcBallCameraControls::onCursorMovement( double x,
                                              double y )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Flip y
    transformScreenCoordinates( x, y, this->_p_window );

    // Update
    this->_dx = x - this->_x;
    this->_dy = y - this->_y;
    this->_x  = x;
    this->_y  = y;

    // Rotate
    if ( this->_isLeftMouseButtonPressed || this->_isRightMouseButtonPressed )
        this->rotate();

    CIE_END_EXCEPTION_TRACING
}


void ArcBallCameraControls::onVerticalScroll( double offset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( offset > 0.0 )
        this->zoomIn();
    else if ( offset < 0.0 )
        this->zoomOut();

    CIE_END_EXCEPTION_TRACING
}


ArcBallCameraControls::configuration_map_type
ArcBallCameraControls::makeConfigurationMap( ArcBallCameraControls::configuration_contents& r_configContents )
{
    CIE_BEGIN_EXCEPTION_TRACING

    ArcBallCameraControls::configuration_map_type configMap;

    auto findControl = [&r_configContents, this]( const std::string& r_controlName )
    -> typename ArcBallCameraControls::configuration_contents::iterator
    {
        auto it = std::find_if(
            r_configContents.begin(),
            r_configContents.end(),
            [&r_controlName]( typename ArcBallCameraControls::configuration_contents::value_type& r_pair )
            {
                if ( r_controlName == r_pair.first )
                    return true;
                return false;
            }
        );

        if ( it == r_configContents.end() )
            this->_p_camera->log( "Control '" + r_controlName + "' not found in the configuration file",
                                  LOG_TYPE_ERROR );

        return it;
    };

    auto it_pair = findControl( "moveForward" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::moveForward, this )
    );

    it_pair = findControl( "moveBackward" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::moveBackward, this )
    );

    it_pair = findControl( "rotateLeft" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::rotateLeft, this )
    );

    it_pair = findControl( "rotateRight" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::rotateRight, this )
    );

    it_pair = findControl( "rotateUp" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::rotateUp, this )
    );

    it_pair = findControl( "rotateDown" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::rotateDown, this )
    );

    it_pair = findControl( "increaseZoomScale" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::increaseZoomScale, this )
    );

    it_pair = findControl( "decreaseZoomScale" );
    configMap.emplace(
        it_pair->second,
        std::bind( &ArcBallCameraControls::decreaseZoomScale, this )
    );

    return configMap;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl
