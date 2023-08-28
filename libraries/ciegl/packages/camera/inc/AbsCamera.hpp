#ifndef CIE_GL_ABS_CAMERA_HPP
#define CIE_GL_ABS_CAMERA_HPP

// --- External Includes ---
#include <glm/glm.hpp>

// --- Utility Includes ---
#include "packages/logging/inc/Loggee.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/RigidBody.hpp"

// --- STL Includes ---
#include <utility>
#include <memory>


namespace cie::gl {


/**
 * Rigid body with view and transformation matrices
 *
 * The projection matrix and its functionalities must
 * be implemented by derived classes
 */
class AbsCamera :
    public utils::Loggee,
    public RigidBody
{
public:
    using internal_matrix_type = glm::mat4;
    using internal_vector_type = glm::vec4;

    using external_vector_type = glm::dvec3;

    CIE_DEFINE_CLASS_POINTERS(AbsCamera)

public:
    AbsCamera( utils::Logger& r_logger                   = utils::LoggerSingleton::get(),
               const std::string& r_name                 = "AbsCamera",
               const AbsCamera::vector_type& r_position  = {0.0, 0.0, 1.0},
               const AbsCamera::vector_type& r_direction = {0.0, 0.0, -1.0},
               const AbsCamera::vector_type& r_up        = {0.0, 1.0, 0.0} );

    AbsCamera( const AbsCamera& r_rhs ) = default;

    AbsCamera( AbsCamera&& r_rhs ) = default;

    AbsCamera& operator=( const AbsCamera& r_rhs ) = delete;

    /**
     * Update matrices
     * By default, transformation matrices are kept
     * up-to-date internally, but derived classes might
     * need extra updating.
     */
    virtual void update();

    /**
     * Set position
     * Trigger view and transformation matrix updates
     */
    virtual void setPosition( const AbsCamera::vector_type& r_position ) override;

    /**
     * Set position and orientation
     * Trigger view and transformation matrix updates
     */
    virtual void setPose( const AbsCamera::vector_type& r_position,
                          const AbsCamera::vector_type& r_direction,
                          const AbsCamera::vector_type& r_up ) override;

    /**
     * Rotate around axis going through the origin
     * Trigger view and transformation matrix updates
     */
    virtual void rotate( double radians,
                         const AbsCamera::vector_type& r_axis ) override;

    /**
     * Rotate around an axis,
     * Trigger view and transformation matrix updates
     */
    virtual void rotate( double radians,
                         const AbsCamera::vector_type& r_axis,
                         const AbsCamera::vector_type& r_pointOnAxis ) override;

    /**
     * Get up-to-date view matrix
     */
    const internal_matrix_type& viewMatrix() const;

    /**
     * Get up-to-date projection matrix
     */
    virtual const internal_matrix_type& projectionMatrix() const = 0;

    /**
     * Get up-to-date transformation matrix
     */
    const internal_matrix_type& transformationMatrix() const;

    /**
     * Get near and far clipping planes
     */
    virtual const std::pair<double,double>& clippingPlanes() const = 0;

    /**
     * Set near and far clipping planes
     * Trigger projection and transformation matrix updates
     */
    virtual void setClippingPlanes( double nearClippingPlane,
                                    double farClippingPlane ) = 0;

    /**
     * Camera type specific zooming
     * Trigger projection and transformation matrix updates
     */
    virtual void zoom( double scale ) = 0;

    virtual const double aspectRatio() const = 0;

    /**
     * Set camera aspect ratio
     * Trigger projection and transformation matrix updates
     */
    virtual void setAspectRatio( double aspectRatio ) = 0;

    /**
     * Field of view in radians (vertical FoV)
     */
    virtual const double fieldOfView() const = 0;

    /**
     * Field of view in radians (vertical FoV)
     * Trigger projection and transformation matrix updates
     */
    virtual void setFieldOfView( double radians ) = 0;

protected:
    /**
     * Rotate around axis going through the camera
     * Trigger view and transformation matrix updates
     */
    virtual void rotateInPlace( double radians,
                                const AbsCamera::vector_type& r_axis ) override;

    /**
     * Recompute view matrix
     */
    void updateViewMatrix();

    /**
     * Recompute projection matrix
     */
    virtual void updateProjectionMatrix() = 0;

    /**
     * Recompute transformation matrix
     * Expects view and projection matrices to be up-to-date
     */
    void updateTransformationMatrix();

protected:
    internal_matrix_type _viewMatrix;
    internal_matrix_type _transformationMatrix;
};


using CameraPtr      = std::shared_ptr<AbsCamera>;
using CameraConstPtr = std::shared_ptr<const AbsCamera>;


} // namespace cie::gl

#endif