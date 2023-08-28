#define _USE_MATH_DEFINES

// --- External Includes ---
#include <glm/gtc/matrix_transform.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/OrthographicProjection.hpp"

// --- STL Includes ---
#include <math.h>


namespace cie::gl {


void OrthographicProjection::zoom_impl( double scale )
{
    CIE_DEBUG_CHECK(
        scale > 0,
        "Invalid zoom scale: " + std::to_string(scale)
    )

    this->_fieldOfView *= 1.0 / scale;
}


void OrthographicProjection::updateProjectionMatrix_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto tmpSizes   = this->sizes();
    tmpSizes.first  /= 2.0;
    tmpSizes.second /= 2.0;

    this->_projectionMatrix = glm::ortho(
        -tmpSizes.first
        ,tmpSizes.first
        ,-tmpSizes.second
        ,tmpSizes.second
        ,this->_clippingPlanes.first
        ,this->_clippingPlanes.second
    );

    CIE_END_EXCEPTION_TRACING
}


std::pair<double,double> OrthographicProjection::sizes() const
{
    std::pair<double,double> sizes;

    sizes.second = 2.0 * this->_clippingPlanes.second
                   * std::tan( this->_fieldOfView * M_PI / 360.0 );
    sizes.first  = this->_aspectRatio * sizes.second;

    return sizes;
}


} // namespace cie::gl