#ifndef CIE_GL_ORTHOGRAPHIC_PROJECTION_HPP
#define CIE_GL_ORTHOGRAPHIC_PROJECTION_HPP

// --- Internal Includes ---
#include "packages/camera/inc/ProjectionPolicy.hpp"

// --- STL Includes ---
#include <utility>


namespace cie::gl {


class OrthographicProjection : public ProjectionPolicy
{
protected:
    void zoom_impl( double scale ) override;
    void updateProjectionMatrix_impl() override;

private:
    /**
     * Projection width and heigh computed from the
     * near clipping plane distance and the vertical FoV
     * 
     * (Since the FoV is undefined for orthographic projections,
     *  it's used to keep track of the projection width and height)
     */
    std::pair<double,double> sizes() const;
};


} // namespace cie::gl


#endif