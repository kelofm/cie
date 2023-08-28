#ifndef CIE_GL_PERSPECTIVE_CAMERA_HPP
#define CIE_GL_PERSPECTIVE_CAMERA_HPP

// --- GL Includes ---
#include "packages/camera/inc/Camera.hpp"
#include "packages/camera/inc/PerspectiveProjection.hpp"


namespace cie::gl {


using PerspectiveCamera = Camera<PerspectiveProjection>;


} // namespace cie::gl

#endif