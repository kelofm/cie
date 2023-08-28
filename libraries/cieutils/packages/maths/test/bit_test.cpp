// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/bit.hpp"


namespace cie::utils {


CIE_TEST_CASE( "flipBit", "[maths]" )
{
    CIE_TEST_CASE_INIT( "flipBit" )

    Size test = 0;

    CIE_TEST_CHECK_NOTHROW( test = flipBit(test, 3) );
    CIE_TEST_CHECK( test == 8 );

    CIE_TEST_CHECK_NOTHROW( test = flipBit(test, 10) );
    CIE_TEST_CHECK( test == 8 + 1024 );

    CIE_TEST_CHECK_NOTHROW( test = flipBit(test, 3) );
    CIE_TEST_CHECK( test == 1024 );
}


CIE_TEST_CASE("getMSBIndex", "[maths]")
{
    CIE_TEST_CASE_INIT("getMSBIndex")

    CIE_TEST_CHECK(getMSBIndex(static_cast<unsigned int>(0)) == 0);
    CIE_TEST_CHECK(getMSBIndex(static_cast<unsigned int>(1)) == 0);
    CIE_TEST_CHECK(getMSBIndex(static_cast<unsigned int>(2)) == 1);
    CIE_TEST_CHECK(getMSBIndex(static_cast<unsigned int>(3)) == 1);
    CIE_TEST_CHECK(getMSBIndex(static_cast<unsigned int>(127)) == 6);
    CIE_TEST_CHECK(getMSBIndex(~static_cast<unsigned int>(0)) == bitSize<unsigned int> - 1);
}


CIE_TEST_CASE("getNumberOfTrailingZeros", "[maths]")
{
    CIE_TEST_CASE_INIT("getNumberOfTrailingZeros")

    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(0)) == bitSize<unsigned int>);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(1)) == 0);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(2)) == 1);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(4)) == 2);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(8)) == 3);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(static_cast<unsigned int>(16)) == 4);
    CIE_TEST_CHECK(getNumberOfTrailingZeros(~static_cast<unsigned int>(127)) == 7);
}


} // namespace cie::utils