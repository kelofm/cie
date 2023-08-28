// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"


namespace cie {


CIE_TEST_CASE( "DynamicArray", "[stl_extension]" )
{
    CIE_TEST_CASE_INIT( "DynamicArray" )

    using Vector = DynamicArray<Size>;

    CIE_TEST_CHECK_NOTHROW( Vector() );
    Vector vector;
    CIE_TEST_REQUIRE( vector.size() == 0 );
    CIE_TEST_CHECK( vector.begin() == vector.end() );

    CIE_TEST_CHECK_NOTHROW( Vector(1) );
    vector = Vector(1);
    CIE_TEST_REQUIRE( vector.size() == 1 );
    CIE_TEST_CHECK( vector.begin() + 1 == vector.end() );
    CIE_TEST_CHECK_NOTHROW( vector.at(0) = 1 );
    CIE_TEST_CHECK( vector[0] == 1 );

    CIE_TEST_CHECK_NOTHROW( Vector(10) );
    vector = Vector( 10 );
    CIE_TEST_REQUIRE( vector.size() == 10 );

    for ( Size i=0; i<vector.size(); ++i )
        vector[i] = i + 1;

    for ( Size i=0; i<vector.size(); ++i )
        CIE_TEST_CHECK( vector[i] == i + 1 );

    CIE_TEST_CHECK_NOTHROW( vector.push_back(11) );
    CIE_TEST_REQUIRE( vector.size() == 11 );
    CIE_TEST_CHECK( vector[10] == 11 );

    CIE_TEST_CHECK_NOTHROW( vector.emplace_back(12) );
    CIE_TEST_REQUIRE( vector.size() == 12 );
    CIE_TEST_CHECK( vector[11] == 12 );

    CIE_TEST_CHECK_NOTHROW( Vector(10, 1) );
}


} // namespace cie