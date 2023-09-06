#define _USE_MATH_DEFINES

// --- External Includes ---
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/control/inc/PanZoomCameraControls.hpp"
#include "packages/utility/inc/glm_overloads.hpp"
#include "packages/control/inc/utilities.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <fstream>
#include <algorithm>
#include <functional>
#include <math.h>


namespace cie::gl {


PanZoomCameraControls::PanZoomCameraControls( bool verbose ) :
    MappedCameraControls(),
    _movementScale( 1.0 ),
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
    std::filesystem::path confiFilePath = getDataPath() / "ciegl/configuration/controls/PanZoomCameraControls.json";
    std::ifstream configFile( confiFilePath );

    if ( !configFile.is_open() )
        CIE_THROW( Exception, "Failed to open control config file: " + confiFilePath.string() );

    this->configure( configFile );

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::setMovementScale( double scale )
{
    if ( scale < 0.0 )
        this->_p_camera->log( "Movement scale must be non-negative (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );

    this->_movementScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set movement scale to " + std::to_string( scale );
}


void PanZoomCameraControls::setZoomScale( double scale )
{
    if ( scale <= 0.0 )
        this->_p_camera->log( "Zoom scale must be positive (" + std::to_string(scale) + ")",
                              LOG_TYPE_ERROR );
    this->_zoomScale = scale;

    if ( this->_verbose )
        *this->_p_camera << "Set zoom scale to " + std::to_string( scale );
}


void PanZoomCameraControls::moveHorizontally( double amount )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto unitVectorRight = glm::cross(
        this->_p_camera->direction(),
        this->_p_camera->up()
    );

    this->_p_camera->translate(
        this->_movementScale * amount * unitVectorRight
    );

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::moveVertically( double amount )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->translate(
        this->_movementScale * amount * this->_p_camera->up()
    );

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::moveLeft()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->moveHorizontally( -1.0 );

    if ( this->_verbose )
        *this->_p_camera << "Move left to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::moveRight()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->moveHorizontally( 1.0 );

    if ( this->_verbose )
        *this->_p_camera << "Move right to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::moveUp()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->moveVertically( 1.0 );

    if ( this->_verbose )
        *this->_p_camera << "Move up to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::moveDown()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->moveVertically( -1.0 );

    if ( this->_verbose )
        *this->_p_camera << "Move down to " << this->_p_camera->position();

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::zoomIn()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom in";

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::zoomOut()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_camera->zoom( 1.0 / this->_zoomScale );

    if ( this->_verbose )
        *this->_p_camera << "Zoom out";

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::onMouseButtonPress(KeyEnum button,
                                               [[maybe_unused]] KeyEnum modifiers)
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = true;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = true;
}


void PanZoomCameraControls::onMouseButtonRelease(KeyEnum button,
                                                 [[maybe_unused]] KeyEnum modifiers)
{
    if ( button == GLFW_MOUSE_BUTTON_1 )
        this->_isLeftMouseButtonPressed = false;
    else if ( button == GLFW_MOUSE_BUTTON_2 )
        this->_isRightMouseButtonPressed = false;
}


void PanZoomCameraControls::onCursorMovement( double x,
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

    // Move
    const auto& windowSize = this->getWindow()->getSize();

    if ( this->_isLeftMouseButtonPressed || this->_isRightMouseButtonPressed )
    {
        // TODO: screen-grab effect
        this->moveHorizontally( this->_dx / windowSize.first / this->_movementScale );
        this->moveVertically( this->_dy / windowSize.second / this->_movementScale );
    }

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::onHorizontalScroll( double offset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( offset > 0.0 )
        this->moveRight();
    else if ( offset < 0.0 )
        this->moveLeft();

    CIE_END_EXCEPTION_TRACING
}


void PanZoomCameraControls::onVerticalScroll( double offset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( offset > 0.0 )
        this->zoomIn();
    else if ( offset < 0.0 )
        this->zoomOut();

    CIE_END_EXCEPTION_TRACING
}


PanZoomCameraControls::configuration_map_type
PanZoomCameraControls::makeConfigurationMap( PanZoomCameraControls::configuration_contents& r_configContents )
{
    CIE_BEGIN_EXCEPTION_TRACING

    PanZoomCameraControls::configuration_map_type configMap;

    auto findControl = [&r_configContents, this]( const std::string& r_controlName )
    -> typename PanZoomCameraControls::configuration_contents::iterator
    {
        auto it = std::find_if(
            r_configContents.begin(),
            r_configContents.end(),
            [&r_controlName]( typename PanZoomCameraControls::configuration_contents::value_type& r_pair )
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

    auto it_pair = findControl( "moveLeft" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::moveLeft, this )
    );

    it_pair = findControl( "moveRight" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::moveRight, this )
    );

    it_pair = findControl( "moveUp" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::moveUp, this )
    );

    it_pair = findControl( "moveDown" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::moveDown, this )
    );

    it_pair = findControl( "zoomIn" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::zoomIn, this )
    );

    it_pair = findControl( "zoomOut" );
    configMap.emplace(
        it_pair->second,
        std::bind( &PanZoomCameraControls::zoomOut, this )
    );

    return configMap;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl
