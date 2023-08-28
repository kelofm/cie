#ifndef CIE_GL_CAMERA_UTILITIES_HPP
#define CIE_GL_CAMERA_UTILITIES_HPP

// --- GEO Includes ---
#include "packages/primitives/inc/Box.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/AbsCamera.hpp"


namespace cie::gl {


template <concepts::Numeric CoordinateType>
void fitViewToBox(AbsCamera& r_camera,
                  const geo::Box<3,CoordinateType>& r_box);


} // namespace cie::gl

#include "packages/camera/impl/camera_utilities_impl.hpp"

#endif