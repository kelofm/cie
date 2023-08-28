// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"


namespace cie::linalg {


CIE_TEST_CASE( "DynamicEigenMatrix", "[matrix]" )
{
    CIE_TEST_CASE_INIT( "DynamicEigenMatrix" )

    DynamicEigenMatrix<int> matrix0, matrix1;

    CIE_TEST_CHECK_NOTHROW( matrix1.resize(2, 3) );

    CIE_TEST_REQUIRE( matrix1.rowSize() == 2 );
    CIE_TEST_REQUIRE( matrix1.columnSize() == 3 );

    CIE_TEST_CHECK_NOTHROW( matrix1(0, 0) = 0 );
    CIE_TEST_CHECK_NOTHROW( matrix1(0, 1) = 1 );
    CIE_TEST_CHECK_NOTHROW( matrix1(0, 2) = 2 );
    CIE_TEST_CHECK_NOTHROW( matrix1(1, 0) = 3 );
    CIE_TEST_CHECK_NOTHROW( matrix1(1, 1) = 4 );
    CIE_TEST_CHECK_NOTHROW( matrix1(1, 2) = 5 );

    CIE_TEST_CHECK( matrix1(0, 0) == 0 );
    CIE_TEST_CHECK( matrix1(0, 1) == 1 );
    CIE_TEST_CHECK( matrix1(0, 2) == 2 );
    CIE_TEST_CHECK( matrix1(1, 0) == 3 );
    CIE_TEST_CHECK( matrix1(1, 1) == 4 );
    CIE_TEST_CHECK( matrix1(1, 2) == 5 );

    CIE_TEST_CHECK_NOTHROW( matrix0 = matrix1 );

    CIE_TEST_REQUIRE( matrix0.rowSize() == 2 );
    CIE_TEST_REQUIRE( matrix0.columnSize() == 3 );

    CIE_TEST_CHECK( matrix0(0, 0) == 0 );
    CIE_TEST_CHECK( matrix0(0, 1) == 1 );
    CIE_TEST_CHECK( matrix0(0, 2) == 2 );
    CIE_TEST_CHECK( matrix0(1, 0) == 3 );
    CIE_TEST_CHECK( matrix0(1, 1) == 4 );
    CIE_TEST_CHECK( matrix0(1, 2) == 5 );

    CIE_TEST_CHECK_NOTHROW( matrix0.resize(3, 2) );
    CIE_TEST_CHECK_NOTHROW( matrix0(0, 0) = 0 );
    CIE_TEST_CHECK_NOTHROW( matrix0(0, 1) = 1 );
    CIE_TEST_CHECK_NOTHROW( matrix0(1, 0) = 2 );
    CIE_TEST_CHECK_NOTHROW( matrix0(1, 1) = 3 );
    CIE_TEST_CHECK_NOTHROW( matrix0(2, 0) = 4 );
    CIE_TEST_CHECK_NOTHROW( matrix0(2, 1) = 5 );
}


} // namespace cie::linalg