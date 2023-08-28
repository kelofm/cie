// --- Internal Includes ---
#include "packages/camera/inc/ProjectionPolicy.hpp"


namespace cie::gl {


ProjectionPolicy::ProjectionPolicy() :
    _projectionMatrix( 1.0 ),
    _clippingPlanes( 1.0, 2.0 ),
    _aspectRatio( 1.0 ),
    _fieldOfView( 60.0 )
{
}


} // namespace cie::gl