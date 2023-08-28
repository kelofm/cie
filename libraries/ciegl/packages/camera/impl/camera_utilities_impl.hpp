#ifndef CIE_GL_CAMERA_UTILITIES_IMPL_HPP
#define CIE_GL_CAMERA_UTILITIES_IMPL_HPP

// --- External Includes ---
#include "glm/glm.hpp"
#include "glm/gtx/rotate_vector.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <concepts::Numeric CoordinateType>
void fitViewToBox(AbsCamera& r_camera,
                  const geo::Box<3,CoordinateType>& r_box)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Get camera direction
    AbsCamera::vector_type boxCenter {
        r_box.base()[0] + r_box.lengths()[0] / 2.0,
        r_box.base()[1] + r_box.lengths()[1] / 2.0,
        r_box.base()[2] + r_box.lengths()[2] / 2.0,
    };

    AbsCamera::vector_type direction = glm::normalize(
        boxCenter - r_camera.position()
    );

    // Get camera up
    AbsCamera::vector_type up;
    if ( std::abs(direction[0]) < 1e-10
         && std::abs(direction[1]) < 1e-10
         && std::abs(direction[2]*direction[2] - 1.0) < 1e-10 )
        up = { 0.0, 1.0, 0.0 };
    else
        up = { 0.0, 0.0, 1.0 };

    up = glm::normalize(glm::cross(up,
                                   direction));

    up = glm::rotate(direction,
                     M_PI / 2.0,
                     up);

    r_camera.setPose(r_camera.position(),
                     direction,
                     up);

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif