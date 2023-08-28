#define _USE_MATH_DEFINES

// --- External Includes ---
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/control/inc/FlyCameraControls.hpp"
#include "packages/utility/inc/glm_overloads.hpp"
#include "packages/control/inc/utilities.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <fstream>
#include <algorithm>
#include <functional>
#include <math.h>


namespace cie::gl {


FlyCameraControls::FlyCameraControls( bool verbose ) :
    MappedCameraControls(),
    _movementScale( 1.0 ),
    _rotationScale( 0.5 * M_PI / 180.0 ),
    _zoomScale( 1.1 ),
    _x( 0.0 ),
    _y( 0.0 ),
    _dx( 0.0 ),
    _dy( 0.0 ),
    _isLeftMouseButtonPressed( false ),
    _isRightMouseButtonPressed( false ),
    _verbose( verbose )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Set default controls
    std::filesystem::path confiFilePath = getDataPath() / "ciegl/configuration/controls/FlyCameraControls.json";
    std::ifstream configFile( confiFilePath );

    if ( !configFile.is_open() )
        CIE_THROW( Exception, "Failed to open control config file: " + confiFilePath.string() );

    this->configure( configFile );

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::setMovementScale( double scale )
{
    if ( scale < 0.0 )
        this->_p_camera->log( "Movement scale must be non-negative (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );

    this->_movementScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set movement scale to " + std::to_string( scale );
}


void FlyCameraControls::setRotationScale( double scale )
{
    if ( scale < 0.0 )
        this->_p_camera->log( "Rotation scale must be non-negative (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );

    this->_rotationScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set rotation scale to " + std::to_string( scale );
}


void FlyCameraControls::setZoomScale( double scale )
{
    if ( scale <= 0.0 )
        this->_p_camera->log( "Zoom scale must be positive (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );
    this->_zoomScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set zoom scale to " + std::to_string( scale );
}


void FlyCameraControls::moveForward()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        this->_movementScale * this->_p_camera->direction()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move forward to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::moveBackward()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        -this->_movementScale * this->_p_camera->direction()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move backward to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::strafeLeft()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto unitVectorRight = glm::cross(
        this->_p_camera->direction(),
        this->_p_camera->up()
    );

    this->_p_camera->translate(
        -this->_movementScale * unitVectorRight
    );

    if ( this->_verbose )
        *this->_p_camera << "Strafe left to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::strafeRight()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto unitVectorRight = glm::cross(
        this->_p_camera->direction(),
        this->_p_camera->up()
    );

    this->_p_camera->translate(
        this->_movementScale * unitVectorRight
    );

    if ( this->_verbose )
        *this->_p_camera << "Strafe right to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::moveUp()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        this->_movementScale * this->_p_camera->up()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move up to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::moveDown()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        -this->_movementScale * this->_p_camera->up()
    );

    if ( this->_verbose )
        *this->_p_camera << "Move down to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::rotate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( this->_dx!=0.0 || this->_dy!=0.0 )
    {
        auto unitVectorRight = glm::cross(
            this->_p_camera->direction(),
            this->_p_camera->up()
        );

        auto axis = this->_dy * unitVectorRight
                    - this->_dx * this->_p_camera->up();

        this->_p_camera->rotate( this->_rotationScale,
                                 axis,
                                 this->_p_camera->position() );
    }

    if ( this->_verbose )
        *this->_p_camera << "Rotate view";

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::rollCounterClockwise()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->rotateRoll( -this->_rotationScale );

    if ( this->_verbose )
        *this->_p_camera << "Roll counter clockwise";

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::rollClockwise()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->rotateRoll( this->_rotationScale );

    if ( this->_verbose )
        *this->_p_camera << "Roll clockwise";

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::zoomIn()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom in";

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::zoomOut()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( 1.0 / this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom out";

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::onMouseButtonPress( KeyEnum button,
                                            KeyEnum modifiers )
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = true;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = true;
}


void FlyCameraControls::onMouseButtonRelease( KeyEnum button,
                                              KeyEnum modifiers )
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = false;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = false;
}


void FlyCameraControls::onCursorMovement( double x,
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


void FlyCameraControls::onHorizontalScroll( double offset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( offset > 0.0 )
        this->strafeRight();
    else if ( offset < 0.0 )
        this->strafeLeft();

    CIE_END_EXCEPTION_TRACING
}


void FlyCameraControls::onVerticalScroll( double offset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( offset > 0.0 )
        this->zoomIn();
    else if ( offset < 0.0 )
        this->zoomOut();

    CIE_END_EXCEPTION_TRACING
}


FlyCameraControls::configuration_map_type
FlyCameraControls::makeConfigurationMap( FlyCameraControls::configuration_contents& r_configContents )
{
    CIE_BEGIN_EXCEPTION_TRACING

    FlyCameraControls::configuration_map_type configMap;

    auto findControl = [&r_configContents, this]( const std::string& r_controlName )
    -> typename FlyCameraControls::configuration_contents::iterator
    {
        auto it = std::find_if(
            r_configContents.begin(),
            r_configContents.end(),
            [&r_controlName]( typename FlyCameraControls::configuration_contents::value_type& r_pair )
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
        std::bind( &FlyCameraControls::moveForward, this )
    );

    it_pair = findControl( "moveBackward" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::moveBackward, this )
    );

    it_pair = findControl( "strafeLeft" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::strafeLeft, this )
    );

    it_pair = findControl( "strafeRight" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::strafeRight, this )
    );

    it_pair = findControl( "moveUp" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::moveUp, this )
    );

    it_pair = findControl( "moveDown" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::moveDown, this )
    );

    it_pair = findControl( "rollCounterClockwise" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::rollCounterClockwise, this )
    );

    it_pair = findControl( "rollClockwise" );
    configMap.emplace(
        it_pair->second,
        std::bind( &FlyCameraControls::rollClockwise, this )
    );

    return configMap;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl