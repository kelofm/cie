// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/overloads/inc/vectoroperators.hpp"

namespace cie {

CIE_TEST_CASE( "DoubleVector operator overloads", "[overloads]" )
{
    DoubleVector v1({ 0.0, 1.0, 2.0 });
    DoubleVector v2({-5.0, 0.0, 5.0 });
    DoubleVector v3({ 10.0, 15.0, 14.0, 12.0});
    DoubleVector resultDoubleVector;
    double scalar = 2.5;
    double resultScalar;

    // SCALAR ADDITION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 + scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + scalar );
    }

    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = scalar + v1 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + scalar );
    }

    // SCALAR SUBTRACTION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 - scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] - scalar );
    }

    // SCALAR MULTIPLICATION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 * scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] * scalar );
    }

    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = scalar * v1 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] * scalar );
    }

    // SCALAR DIVISION
    #ifdef CIE_ENABLE_DIVISION_BY_ZERO_CHECKS
        CIE_TEST_REQUIRE_THROWS( resultDoubleVector = v1 / 0.0 );
    #endif
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 / scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] / scalar );
    }

    // VECTOR ADDITION
    #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        CIE_TEST_REQUIRE_THROWS( resultDoubleVector = v1 + v3 );
    #endif
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 + v2 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + v2[i] );
    }

    // VECTOR SUBTRACTION
    #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        CIE_TEST_REQUIRE_THROWS( resultDoubleVector = v1 - v3 );
    #endif
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 - v2 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] - v2[i] );
    }

    // DOT PRODUCT
    CIE_TEST_REQUIRE_THROWS( resultScalar = v1 * v3 );
    CIE_TEST_REQUIRE_NOTHROW( resultScalar = v1 * v2 );
    CIE_TEST_CHECK( resultScalar == v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2] );
}


CIE_TEST_CASE( "DoubleArray operator overloads", "[overloads]" )
{
    DoubleArray<3> v1   = { 0.0, 1.0, 2.0 };
    DoubleArray<3> v2   = {-5.0, 0.0, 5.0 };
    DoubleArray<3> resultDoubleVector;
    double scalar = 2.5;
    double resultScalar;

    // SCALAR ADDITION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 + scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + scalar );
    }

    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = scalar + v1 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + scalar );
    }

    // SCALAR SUBTRACTION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 - scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] - scalar );
    }

    // SCALAR MULTIPLICATION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 * scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] * scalar );
    }

    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = scalar * v1 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] * scalar );
    }

    // SCALAR DIVISION
    #ifdef CIE_ENABLE_DIVISION_BY_ZERO_CHECKS
        CIE_TEST_REQUIRE_THROWS( resultDoubleVector = v1 / 0.0 );
    #endif
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 / scalar );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] / scalar );
    }

    // VECTOR ADDITION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 + v2 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] + v2[i] );
    }

    // VECTOR SUBTRACTION
    CIE_TEST_REQUIRE_NOTHROW( resultDoubleVector = v1 - v2 );
    CIE_TEST_REQUIRE( resultDoubleVector.size() == v1.size() );
    for (size_t i=0; i<resultDoubleVector.size(); ++i){
        CIE_TEST_CHECK( resultDoubleVector[i]== v1[i] - v2[i] );
    }

    // DOT PRODUCT
    CIE_TEST_REQUIRE_NOTHROW( resultScalar = v1 * v2 );
    CIE_TEST_CHECK( resultScalar == v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2] );
}

} // namespace cie