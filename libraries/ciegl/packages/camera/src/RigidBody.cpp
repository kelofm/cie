// --- External Includes ---
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/RigidBody.hpp"


namespace cie::gl {


RigidBody::RigidBody(   const RigidBody::vector_type& r_position,
                        const RigidBody::vector_type& r_direction,
                        const RigidBody::vector_type& r_up ) :
    _position( r_position ),
    _direction( r_direction ),
    _up( r_up )
{
}


RigidBody::RigidBody() :
    RigidBody( RigidBody::vector_type( 0.0f, 0.0f, 0.0f ),
               RigidBody::vector_type( 1.0f, 0.0f, 0.0f ),
               RigidBody::vector_type( 0.0f, 0.0f, 1.0f ) )
{
}


const RigidBody::vector_type& RigidBody::position() const
{
    return _position;
}


const RigidBody::vector_type& RigidBody::direction() const
{
    return _direction;
}


const RigidBody::vector_type& RigidBody::up() const
{
    return _up;
}


void RigidBody::setPosition( const RigidBody::vector_type& r_position )
{
    _position = r_position;
}


void RigidBody::setPose( const RigidBody::vector_type& r_position,
                         const RigidBody::vector_type& r_direction,
                         const RigidBody::vector_type& r_up )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check whether direction and up are perpendicular
    CIE_DEBUG_CHECK(
        std::abs(glm::dot( r_direction, r_up )) < 1e-8,
        "Camera direction and up do not define a cartesian system!"
    )

    _position   = r_position;
    _direction  = glm::normalize( r_direction );
    _up         = glm::normalize( r_up );

    CIE_END_EXCEPTION_TRACING
}


void RigidBody::translate( const RigidBody::vector_type& r_translation )
{
    this->setPosition( this->_position + r_translation );
}


void RigidBody::rotate( double radians,
                        const RigidBody::vector_type& r_axis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_DEBUG_CHECK(
        std::abs(r_axis[0]) + std::abs(r_axis[1]) + std::abs(r_axis[2]) > 1e-15,
        "Rotation axis is a nullvector"
    )

    RigidBody::matrix_type rotation = glm::rotate( radians, r_axis );
    _position  = rotation * glm::dvec4( _position, 1.0 );
    _direction = rotation * glm::dvec4( _direction, 1.0 );
    _up        = rotation * glm::dvec4( _up, 1.0 );

    CIE_END_EXCEPTION_TRACING
}


void RigidBody::rotate( double radians,
                        const RigidBody::vector_type& r_axis,
                        const RigidBody::vector_type& r_pointOnAxis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_DEBUG_CHECK(
        std::abs(r_axis[0]) + std::abs(r_axis[1]) + std::abs(r_axis[2]) > 1e-15,
        "Rotation axis is a nullvector"
    )

    // Create transformation tensors
    glm::dmat4 rotation = glm::rotate( radians, r_axis );

    // Compute closest point on the axis from the current position
    //  -> body position in the rotation's local coordinate system
    RigidBody::vector_type closestPointOnAxis = r_pointOnAxis
                                                + (glm::dot( r_pointOnAxis, r_axis) - glm::dot( _position, r_axis ) )
                                                  / glm::dot( r_axis, r_axis )
                                                  * r_axis;

    // Rotate body position in the local system,
    // and add it the the position of the local origin
    //  -> rotation around an arbitrary axis
    _position = closestPointOnAxis
                + glm::dvec3( rotation * glm::dvec4( _position - closestPointOnAxis, 1.0 ) );

    // Local axes need only be rotated
    _direction = rotation * glm::dvec4( _direction, 1.0 );
    _up        = rotation * glm::dvec4( _up, 1.0 );

    CIE_END_EXCEPTION_TRACING
}


void RigidBody::rotateYaw( double radians )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->rotateInPlace( radians,
                         this->getLocalZAxis() );

    CIE_END_EXCEPTION_TRACING
}


void RigidBody::rotatePitch( double radians )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->rotateInPlace( radians,
                         this->getLocalYAxis() );

    CIE_END_EXCEPTION_TRACING
}


void RigidBody::rotateRoll( double radians )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->rotateInPlace( radians,
                         this->getLocalXAxis() );

    CIE_END_EXCEPTION_TRACING
}


RigidBody::vector_type RigidBody::getLocalXAxis() const
{
    return _direction;
}


RigidBody::vector_type RigidBody::getLocalYAxis() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return glm::cross( _up, _direction );

    CIE_END_EXCEPTION_TRACING
}


RigidBody::vector_type RigidBody::getLocalZAxis() const
{
    return _up;
}


void RigidBody::rotateInPlace( double radians,
                               const RigidBody::vector_type& r_axis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_DEBUG_CHECK(
        std::abs(r_axis[0]) + std::abs(r_axis[1]) + std::abs(r_axis[2]) > 1e-15,
        "Rotation axis is a nullvector"
    )

    glm::dmat4 rotation = glm::rotate( radians, r_axis );

    _direction          = rotation * glm::dvec4( _direction, 1.0 );
    _up                 = rotation * glm::dvec4( _up, 1.0 );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl