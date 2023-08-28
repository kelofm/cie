// --- External Includes ---
#include <glm/gtc/matrix_transform.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/AbsCamera.hpp"


namespace cie::gl {


AbsCamera::AbsCamera( utils::Logger& r_logger,
                      const std::string& r_name,
                      const AbsCamera::vector_type& r_position,
                      const AbsCamera::vector_type& r_direction,
                      const AbsCamera::vector_type& r_up ) :
    utils::Loggee( r_logger, r_name ),
    RigidBody( r_position, r_direction, r_up ),
    _viewMatrix( 1.0 ),
    _transformationMatrix( 1.0 )
{
    // Derived class constructors should call setPose!
    // Since this class is constructed before its derived ones,
    // setPose would indirectly call pure virtual functions (updateProjectionMatrix)
}


void AbsCamera::update()
{
}


void AbsCamera::setPosition( const AbsCamera::vector_type& r_position )
{
    CIE_BEGIN_EXCEPTION_TRACING

    RigidBody::setPosition( r_position );

    this->updateViewMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}



void AbsCamera::setPose( const AbsCamera::vector_type& r_position,
                         const AbsCamera::vector_type& r_direction,
                         const AbsCamera::vector_type& r_up )
{
    CIE_BEGIN_EXCEPTION_TRACING

    RigidBody::setPose( r_position,
                        r_direction,
                        r_up );
    
    this->updateViewMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


void AbsCamera::rotate( double radians,
                        const AbsCamera::vector_type& r_axis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    RigidBody::rotate( radians, r_axis );

    this->updateViewMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


void AbsCamera::rotate( double radians,
                        const AbsCamera::vector_type& r_axis,
                        const AbsCamera::vector_type& r_pointOnAxis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    RigidBody::rotate( radians,
                       r_axis,
                       r_pointOnAxis );

    this->updateViewMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


void AbsCamera::rotateInPlace( double radians,
                               const AbsCamera::vector_type& r_axis )
{
    CIE_BEGIN_EXCEPTION_TRACING

    RigidBody::rotateInPlace( radians, r_axis );

    this->updateViewMatrix();
    this->updateTransformationMatrix();

    CIE_END_EXCEPTION_TRACING
}


const AbsCamera::internal_matrix_type& AbsCamera::viewMatrix() const
{
    return this->_viewMatrix;
}


const AbsCamera::internal_matrix_type& AbsCamera::transformationMatrix() const
{
    return this->_transformationMatrix;
}


void AbsCamera::updateViewMatrix()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_viewMatrix = glm::lookAt( this->_position,
                                     this->_position + this->_direction,
                                     this->_up );

    CIE_END_EXCEPTION_TRACING
}


void AbsCamera::updateTransformationMatrix()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_transformationMatrix = this->projectionMatrix() * this->_viewMatrix;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl