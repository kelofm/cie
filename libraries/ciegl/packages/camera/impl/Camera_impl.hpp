#ifndef CIE_GL_CAMERA_IMPL_HPP
#define CIE_GL_CAMERA_IMPL_HPP

#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/ProjectionPolicy.hpp"

// --- STL Includes ---
#include <math.h>


namespace cie::gl {


template <concepts::ProjectionPolicy ProjectionType>
Camera<ProjectionType>::Camera( utils::Logger& r_logger,
                                const std::string& r_name,
                                const typename Camera<ProjectionType>::vector_type& r_position,
                                const typename Camera<ProjectionType>::vector_type& r_direction,
                                const typename Camera<ProjectionType>::vector_type& r_up ) :
    AbsCamera( r_logger,
               r_name,
               r_position,
               r_direction,
               r_up ),
    ProjectionType()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->setPose( r_position,
                   r_direction,
                   r_up );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ProjectionPolicy ProjectionType>
inline const typename AbsCamera::internal_matrix_type&
Camera<ProjectionType>::projectionMatrix() const
{
    return this->_projectionMatrix;
}


template <concepts::ProjectionPolicy ProjectionType>
inline const double
Camera<ProjectionType>::fieldOfView() const
{
    return this->_fieldOfView;
}


template <concepts::ProjectionPolicy ProjectionType>
inline void
Camera<ProjectionType>::setFieldOfView( double radians )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( radians < 0 || M_PI < radians ) [[unlikely]]
        CIE_THROW( Exception, "Invalid field of view: " + std::to_string(radians) );

    this->_fieldOfView = radians;
    this->updateProjectionMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ProjectionPolicy ProjectionType>
inline void
Camera<ProjectionType>::zoom( double scale )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->zoom_impl( scale );
    this->updateProjectionMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ProjectionPolicy ProjectionType>
inline void
Camera<ProjectionType>::updateProjectionMatrix()
{
    this->updateProjectionMatrix_impl();
}


template <concepts::ProjectionPolicy ProjectionType>
inline const std::pair<double,double>&
Camera<ProjectionType>::clippingPlanes() const
{
    return this->_clippingPlanes;
}


template <concepts::ProjectionPolicy ProjectionType>
inline void
Camera<ProjectionType>::setClippingPlanes( double nearClippingPlane,
                                           double farClippingPlane )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( nearClippingPlane < 0 || farClippingPlane < nearClippingPlane ) [[unlikely]]
        CIE_THROW(
            Exception,
            "Invalid clipping planes: " + std::to_string(nearClippingPlane)
                + " " + std::to_string(farClippingPlane)
        )

    this->_clippingPlanes.first  = nearClippingPlane;
    this->_clippingPlanes.second = farClippingPlane;

    this->updateProjectionMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ProjectionPolicy ProjectionType>
inline const double
Camera<ProjectionType>::aspectRatio() const
{
    return this->_aspectRatio;
}


template <concepts::ProjectionPolicy ProjectionType>
inline void
Camera<ProjectionType>::setAspectRatio( double aspectRatio )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_DEBUG_CHECK(
        aspectRatio > 0,
        "Invalid aspect ratio: " + std::to_string(aspectRatio)
    )

    this->_aspectRatio = aspectRatio;
    this->updateProjectionMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif