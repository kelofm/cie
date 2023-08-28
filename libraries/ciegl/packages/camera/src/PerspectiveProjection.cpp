// --- External Includes ---
#include <glm/gtc/matrix_transform.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/PerspectiveProjection.hpp"


namespace cie::gl {


void PerspectiveProjection::zoom_impl( double scale )
{
    CIE_DEBUG_CHECK(
        scale > 0,
        "Invalid zoom scale: " + std::to_string(scale)
    )

    this->_fieldOfView *= 1.0 / scale;
}


void PerspectiveProjection::updateProjectionMatrix_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_projectionMatrix = glm::perspective( this->_fieldOfView,
                                                this->_aspectRatio,
                                                this->_clippingPlanes.first,
                                                this->_clippingPlanes.second );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl