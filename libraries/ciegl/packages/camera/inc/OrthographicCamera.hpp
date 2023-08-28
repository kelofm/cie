#ifndef CIE_GL_ORTHOGRAPHIC_CAMERA_HPP
#define CIE_GL_ORTHOGRAPHIC_CAMERA_HPP

// --- GL Includes ---
#include "packages/camera/inc/Camera.hpp"
#include "packages/camera/inc/OrthographicProjection.hpp"


namespace cie::gl {


using OrthographicCamera = Camera<OrthographicProjection>;


} // namespace cie::gl

#endif