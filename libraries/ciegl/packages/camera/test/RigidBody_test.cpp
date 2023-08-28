#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/RigidBody.hpp"

// --- STL Includes ---
#include <math.h>


namespace cie::gl {


CIE_TEST_CASE( "RigidBody", "[rigidbody]" )
{
    CIE_TEST_CASE_INIT( "RigidBody" )
    
    RigidBodyPtr p_body;

    auto resetBody = [&]() -> void
    {
        p_body->setPose(
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        );
    };

    {
        CIE_TEST_CASE_INIT( "construction" )

        CIE_TEST_REQUIRE_NOTHROW( p_body.reset(
            new RigidBody( { 0.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 0.0, 0.0, 1.0 } )
        ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(1.0) );
    }

    {
        CIE_TEST_CASE_INIT( "hard setting" )

        CIE_TEST_CHECK_NOTHROW( p_body->setPosition( {1.0, 1.0, 1.0} ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(1.0) );

        CIE_TEST_CHECK_NOTHROW( p_body->setPose(
            { 1.0, 2.0, 3.0 },
            { 0.0, 1.0, 0.0 },
            { 1.0, 0.0, 0.0 }
        ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(2.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(3.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(0.0) );
    }

    {
        CIE_TEST_CASE_INIT( "translation" )

        resetBody();
        CIE_TEST_CHECK_NOTHROW( p_body->translate( {1.0, 1.0, 1.0} ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(1.0) );
    }

    {
        CIE_TEST_CASE_INIT( "origin-bound rotation" )

        resetBody();
        p_body->translate( {1.0, 1.0, 1.0} );
        CIE_TEST_CHECK_NOTHROW( p_body->rotate( M_PI / 2.0,
                                       {0.0, 0.0, 1.0} ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(-1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(1.0) );

        CIE_TEST_CHECK_NOTHROW( p_body->rotate( M_PI / 2.0,
                                       {1.0, 0.0, 0.0} ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(-1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(-1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(-1.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(0.0).margin(1e-15) );
    }

    {
        CIE_TEST_CASE_INIT( "general rotation" )

        resetBody();
        p_body->translate( {1.0, 1.0, 1.0} );
        CIE_TEST_CHECK_NOTHROW( p_body->rotate( M_PI / 2.0,
                                       {0.0, 1.0, 0.0},
                                       p_body->position() ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(-1.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(0.0).margin(1e-15) );

        CIE_TEST_CHECK_NOTHROW( p_body->rotate( M_PI / 2.0,
                                       {0.0, 1.0, 0.0},
                                       {1.0, 1.0, 0.0} ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(2.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(0.0).margin(1e-15) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(-1.0) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0).margin(1e-15) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(-1.0) );
    }

    {
        CIE_TEST_CASE_INIT( "euler rotations" )

        resetBody();
        p_body->translate( {1.0, 1.0, 1.0} );
        CIE_TEST_CHECK_NOTHROW( p_body->rotateYaw( M_PI / 2.0 ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(0.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(1.0) );

        CIE_TEST_CHECK_NOTHROW( p_body->rotatePitch( M_PI / 2.0 ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(-1.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(0.0).margin(1e-15) );

        CIE_TEST_CHECK_NOTHROW( p_body->rotateRoll( M_PI / 2.0 ) );

        CIE_TEST_CHECK( p_body->position()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[1] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->position()[2] == Approx(1.0) );

        CIE_TEST_CHECK( p_body->direction()[0] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->direction()[2] == Approx(-1.0) );

        CIE_TEST_CHECK( p_body->up()[0] == Approx(1.0) );
        CIE_TEST_CHECK( p_body->up()[1] == Approx(0.0).margin(1e-15) );
        CIE_TEST_CHECK( p_body->up()[2] == Approx(0.0).margin(1e-15) );
    }

}


} // namespace cie::gl