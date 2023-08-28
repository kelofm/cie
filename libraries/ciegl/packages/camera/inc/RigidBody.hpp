#ifndef CIEGL_RIGID_BODY_HPP
#define CIEGL_RIGID_BODY_HPP

// --- External Includes ---
#include <glm/glm.hpp>

// --- STL Includes ---
#include <memory>


namespace cie::gl {


/**
 * Minimal class for keeping track of
 * the position and orientation of a
 * rigid body.
 */
class RigidBody
{
public:
    using vector_type       = glm::dvec3;
    using quaternion_type   = glm::dquat;
    using matrix_type       = glm::dmat4;

public:
    RigidBody(  const vector_type& r_position,
                const vector_type& r_direction,
                const vector_type& r_up );
    RigidBody();

    RigidBody( const RigidBody& r_rhs ) = default;
    
    RigidBody( RigidBody&& r_rhs ) = default;

    const vector_type& position() const;
    const vector_type& direction() const;
    const vector_type& up() const;

    /// Move body to the specified position
    virtual void setPosition( const vector_type& r_position );

    /// Set position and orientation
    virtual void setPose( const vector_type& r_position,
                          const vector_type& r_direction,
                          const vector_type& r_up );

    /// Translate body by the specified vector
    virtual void translate( const vector_type& r_translation );

    /// Rotate body around an origin-bound axis
    virtual void rotate( double radians,
                         const vector_type& r_axis );

    /// Rotate body around a general axis
    virtual void rotate( double radians,
                         const vector_type& r_axis,
                         const vector_type& r_pointOnAxis );
    
    /// Rotate body around its local z axis
    void rotateYaw( double radians );

    /// Rotate body around its local y axis
    void rotatePitch( double radians );

    /// Rotate body around its local x axis
    void rotateRoll( double radians );

    /// Roll axis
    vector_type getLocalXAxis() const;

    /// Pitch axis
    vector_type getLocalYAxis() const;

    /// Yaw axis
    vector_type getLocalZAxis() const;

protected:
    /// Rotate around an axis that goes through the body
    virtual void rotateInPlace( double radians,
                                const vector_type& r_axis );

protected:
    vector_type     _position;
    vector_type     _direction;
    vector_type     _up;
};


using RigidBodyPtr      = std::shared_ptr<RigidBody>;
using RigidBodyConstPtr = std::shared_ptr<const RigidBody>;


} // namespace cie::gl


#endif