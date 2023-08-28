#ifndef CIE_GL_PROJECTION_POLICY_HPP
#define CIE_GL_PROJECTION_POLICY_HPP

// --- External Includes ---
#include <glm/glm.hpp>

// --- STL Includes ---
#include <concepts>
#include <utility>


namespace cie::gl {


class ProjectionPolicy
{
protected:
    using internal_matrix_type = glm::mat4;

protected:
    ProjectionPolicy();

    ProjectionPolicy( const ProjectionPolicy& r_rhs ) = default;

    ProjectionPolicy( ProjectionPolicy&& r_rhs ) = default;

    virtual void zoom_impl( double scale ) = 0;
    virtual void updateProjectionMatrix_impl() = 0;

protected:
    internal_matrix_type      _projectionMatrix;
    std::pair<double,double>  _clippingPlanes;
    double                    _aspectRatio;
    double                    _fieldOfView;
};


} // namespace cie::gl


namespace cie::concepts {
template <class T>
concept ProjectionPolicy
= std::derived_from<T,cie::gl::ProjectionPolicy>;
} // namespace cie::concepts


#endif