#ifndef CIE_GL_CAMERA_HPP
#define CIE_GL_CAMERA_HPP

// --- Internal Includes ---
#include "packages/camera/inc/AbsCamera.hpp"
#include "packages/camera/inc/ProjectionPolicy.hpp"

// --- STL Includes ---
#include <utility>


namespace cie::gl {


template <concepts::ProjectionPolicy ProjectionType>
class Camera :
    public AbsCamera,
    public ProjectionType
{
public:
    using projection_type = ProjectionType;

public:
    Camera( utils::Logger& r_logger                   = utils::LoggerSingleton::get(),
            const std::string& r_name                 = "Camera",
            const AbsCamera::vector_type& r_position  = {0.0, 0.0, 1.0},
            const AbsCamera::vector_type& r_direction = {0.0, 0.0, -1.0},
            const AbsCamera::vector_type& r_up        = {0.0, 1.0, 0.0} );

    Camera( const Camera<ProjectionType>& r_rhs ) = default;

    Camera( Camera<ProjectionType>&& r_rhs ) = default;

    Camera<ProjectionType>& operator=( const Camera<ProjectionType>& r_rhs ) = delete;

    const typename AbsCamera::internal_matrix_type& projectionMatrix() const override;

    void zoom( double scale ) override;

    const std::pair<double,double>& clippingPlanes() const override;

    void setClippingPlanes( double nearClippingPlane,
                            double farClippingPlane ) override;

    double aspectRatio() const override;

    void setAspectRatio( double aspectRatio ) override;

    double fieldOfView() const override;

    void setFieldOfView( double radians ) override;

protected:
    void updateProjectionMatrix() override;
};


} // namespace cie::gl

#include "packages/camera/impl/Camera_impl.hpp"

#endif