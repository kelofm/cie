// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Incldues ---
#include "packages/matrix/inc/EigenArray.hpp"
#include "packages/overloads/inc/matrix_operators.hpp"


namespace cie::linalg {


CIE_TEST_CASE( "matrix operators", "[overloads]" )
{
    CIE_TEST_CASE_INIT( "matrix operators" )

    {
        CIE_TEST_CASE_INIT( "DynamicEigenMatrix" )

        DynamicEigenMatrix<int> matrix( 2, 3 );
        matrix( 0, 0 ) = 0;
        matrix( 0, 1 ) = 1;
        matrix( 0, 2 ) = 2;
        matrix( 1, 0 ) = 3;
        matrix( 1, 1 ) = 4;
        matrix( 1, 2 ) = 5;

        {
            CIE_TEST_CASE_INIT( "scalar * matrix" )

            auto tmp = matrix;

            CIE_TEST_CHECK_NOTHROW( tmp *= 2 );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 2 );
            CIE_TEST_CHECK( tmp(0, 2) == 4 );
            CIE_TEST_CHECK( tmp(1, 0) == 6 );
            CIE_TEST_CHECK( tmp(1, 1) == 8 );
            CIE_TEST_CHECK( tmp(1, 2) == 10 );

            CIE_TEST_CHECK_NOTHROW( tmp = 2 * tmp );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 4 );
            CIE_TEST_CHECK( tmp(0, 2) == 8 );
            CIE_TEST_CHECK( tmp(1, 0) == 12 );
            CIE_TEST_CHECK( tmp(1, 1) == 16 );
            CIE_TEST_CHECK( tmp(1, 2) == 20 );
        } // scalar * matrix

        {
            CIE_TEST_CASE_INIT( "matrix * vector" )

            EigenVector<int> vector3( 3 );
            vector3[0] = 10;
            vector3[1] = 20;
            vector3[2] = 30;

            CIE_TEST_REQUIRE_NOTHROW( matrix * vector3 );
            auto result2 = matrix * vector3;

            CIE_TEST_REQUIRE( result2.size() == 2 );
            CIE_TEST_CHECK( result2[0] == 80 );
            CIE_TEST_CHECK( result2[1] == 260 );
        } // matrix * vector

        {
            CIE_TEST_CASE_INIT( "matrix + matrix" )

            auto tmp = matrix;

            CIE_TEST_CHECK_NOTHROW( tmp += 2 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 3 );
            CIE_TEST_CHECK( tmp(0, 2) == 6 );
            CIE_TEST_CHECK( tmp(1, 0) == 9 );
            CIE_TEST_CHECK( tmp(1, 1) == 12 );
            CIE_TEST_CHECK( tmp(1, 2) == 15 );

            CIE_TEST_CHECK_NOTHROW( tmp -= 3 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 0 );
            CIE_TEST_CHECK( tmp(0, 2) == 0 );
            CIE_TEST_CHECK( tmp(1, 0) == 0 );
            CIE_TEST_CHECK( tmp(1, 1) == 0 );
            CIE_TEST_CHECK( tmp(1, 2) == 0 );

            CIE_TEST_CHECK_NOTHROW( tmp = tmp + matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 1 );
            CIE_TEST_CHECK( tmp(0, 2) == 2 );
            CIE_TEST_CHECK( tmp(1, 0) == 3 );
            CIE_TEST_CHECK( tmp(1, 1) == 4 );
            CIE_TEST_CHECK( tmp(1, 2) == 5 );

            CIE_TEST_CHECK_NOTHROW( tmp = tmp - 2 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == -1 );
            CIE_TEST_CHECK( tmp(0, 2) == -2 );
            CIE_TEST_CHECK( tmp(1, 0) == -3 );
            CIE_TEST_CHECK( tmp(1, 1) == -4 );
            CIE_TEST_CHECK( tmp(1, 2) == -5 );
        } // matrix + matrix

        {
            DynamicEigenMatrix<int> lhs( 3, 2 );
            lhs( 0, 0 ) = 0;
            lhs( 0, 1 ) = 1;
            lhs( 1, 0 ) = 2;
            lhs( 1, 1 ) = 3;
            lhs( 2, 0 ) = 4;
            lhs( 2, 1 ) = 5;
            auto tmp = lhs;

            CIE_TEST_CHECK_NOTHROW( tmp = lhs * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 3 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 3 );
            CIE_TEST_CHECK( tmp(0, 1) == 4 );
            CIE_TEST_CHECK( tmp(0, 2) == 5 );
            CIE_TEST_CHECK( tmp(1, 0) == 9 );
            CIE_TEST_CHECK( tmp(1, 1) == 14 );
            CIE_TEST_CHECK( tmp(1, 2) == 19 );
            CIE_TEST_CHECK( tmp(2, 0) == 15 );
            CIE_TEST_CHECK( tmp(2, 1) == 24 );
            CIE_TEST_CHECK( tmp(2, 2) == 33 );
        }
    } // DynamicEigenMatrix

    {
        CIE_TEST_CASE_INIT( "DynamicEigenMatrix" )

        StaticEigenMatrix<int, 2, 3> matrix;
        matrix( 0, 0 ) = 0;
        matrix( 0, 1 ) = 1;
        matrix( 0, 2 ) = 2;
        matrix( 1, 0 ) = 3;
        matrix( 1, 1 ) = 4;
        matrix( 1, 2 ) = 5;

        {
            CIE_TEST_CASE_INIT( "scalar * matrix" )

            auto tmp = matrix;

            CIE_TEST_CHECK_NOTHROW( tmp *= 2 );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 2 );
            CIE_TEST_CHECK( tmp(0, 2) == 4 );
            CIE_TEST_CHECK( tmp(1, 0) == 6 );
            CIE_TEST_CHECK( tmp(1, 1) == 8 );
            CIE_TEST_CHECK( tmp(1, 2) == 10 );

            CIE_TEST_CHECK_NOTHROW( tmp = 2 * tmp );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 4 );
            CIE_TEST_CHECK( tmp(0, 2) == 8 );
            CIE_TEST_CHECK( tmp(1, 0) == 12 );
            CIE_TEST_CHECK( tmp(1, 1) == 16 );
            CIE_TEST_CHECK( tmp(1, 2) == 20 );
        } // scalar * matrix

        {
            CIE_TEST_CASE_INIT( "matrix * vector" )

            EigenArray<int,3> vector3;
            vector3[0] = 10;
            vector3[1] = 20;
            vector3[2] = 30;

            CIE_TEST_REQUIRE_NOTHROW( matrix * vector3 );
            auto result2 = matrix * vector3;

            CIE_TEST_REQUIRE( result2.size() == 2 );
            CIE_TEST_CHECK( result2[0] == 80 );
            CIE_TEST_CHECK( result2[1] == 260 );
        } // matrix * vector

        {
            CIE_TEST_CASE_INIT( "matrix + matrix" )

            auto tmp = matrix;

            CIE_TEST_CHECK_NOTHROW( tmp += 2 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 3 );
            CIE_TEST_CHECK( tmp(0, 2) == 6 );
            CIE_TEST_CHECK( tmp(1, 0) == 9 );
            CIE_TEST_CHECK( tmp(1, 1) == 12 );
            CIE_TEST_CHECK( tmp(1, 2) == 15 );

            CIE_TEST_CHECK_NOTHROW( tmp -= 3 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 0 );
            CIE_TEST_CHECK( tmp(0, 2) == 0 );
            CIE_TEST_CHECK( tmp(1, 0) == 0 );
            CIE_TEST_CHECK( tmp(1, 1) == 0 );
            CIE_TEST_CHECK( tmp(1, 2) == 0 );

            CIE_TEST_CHECK_NOTHROW( tmp = tmp + matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == 1 );
            CIE_TEST_CHECK( tmp(0, 2) == 2 );
            CIE_TEST_CHECK( tmp(1, 0) == 3 );
            CIE_TEST_CHECK( tmp(1, 1) == 4 );
            CIE_TEST_CHECK( tmp(1, 2) == 5 );

            CIE_TEST_CHECK_NOTHROW( tmp = tmp - 2 * matrix );
            CIE_TEST_REQUIRE( tmp.rowSize() == 2 );
            CIE_TEST_REQUIRE( tmp.columnSize() == 3 );
            CIE_TEST_CHECK( tmp(0, 0) == 0 );
            CIE_TEST_CHECK( tmp(0, 1) == -1 );
            CIE_TEST_CHECK( tmp(0, 2) == -2 );
            CIE_TEST_CHECK( tmp(1, 0) == -3 );
            CIE_TEST_CHECK( tmp(1, 1) == -4 );
            CIE_TEST_CHECK( tmp(1, 2) == -5 );
        } // matrix + matrix

        {
            StaticEigenMatrix<int, 3, 2> lhs;
            lhs( 0, 0 ) = 0;
            lhs( 0, 1 ) = 1;
            lhs( 1, 0 ) = 2;
            lhs( 1, 1 ) = 3;
            lhs( 2, 0 ) = 4;
            lhs( 2, 1 ) = 5;

            StaticEigenMatrix<int, 3, 3> result;

            CIE_TEST_CHECK_NOTHROW( result = lhs * matrix );
            CIE_TEST_REQUIRE( result.rowSize() == 3 );
            CIE_TEST_REQUIRE( result.columnSize() == 3 );
            CIE_TEST_CHECK( result(0, 0) == 3 );
            CIE_TEST_CHECK( result(0, 1) == 4 );
            CIE_TEST_CHECK( result(0, 2) == 5 );
            CIE_TEST_CHECK( result(1, 0) == 9 );
            CIE_TEST_CHECK( result(1, 1) == 14 );
            CIE_TEST_CHECK( result(1, 2) == 19 );
            CIE_TEST_CHECK( result(2, 0) == 15 );
            CIE_TEST_CHECK( result(2, 1) == 24 );
            CIE_TEST_CHECK( result(2, 2) == 33 );
        }
    } // StaticEigenMatrix
}


} // namespace cie::linalg