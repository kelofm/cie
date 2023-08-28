// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/base.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"


namespace cie::utils {


CIE_TEST_CASE( "base", "[maths]" )
{
    CIE_TEST_CASE_INIT( "base" )

    const Size                      NDigits = 4;
    typedef Size                    NT;
    typedef StaticArray<NT,NDigits>  Array;

    NT base;

    {
        CIE_TEST_CASE_INIT( "baseN in-place" )

        NT value    = 13;
        base        = 2;
        Array check = { 1, 0, 1, 1 };
        Array test;

        CIE_TEST_REQUIRE_NOTHROW( baseN( value, base, test ) );
        for (Size i=0; i<check.size(); ++i)
            CIE_TEST_CHECK( test[i] == check[i] );


        base        = 5;
        check       = { 3, 2, 0, 0 };
        CIE_TEST_REQUIRE_NOTHROW( baseN( value, base, test ) );
        for (Size i=0; i<check.size(); ++i)
            CIE_TEST_CHECK( test[i] == check[i] );
    }

    {
        CIE_TEST_CASE_INIT( "baseN" )

        NT value    = 13;
        base        = 2;
        Array check = { 1, 0, 1, 1 };
        Array test;

        CIE_TEST_REQUIRE_NOTHROW( test = baseN<NDigits>( value, base ) );
        for (Size i=0; i<check.size(); ++i)
            CIE_TEST_CHECK( test[i] == check[i] );

        base        = 5;
        check       = { 3, 2, 0, 0 };
        CIE_TEST_REQUIRE_NOTHROW( test = baseN<NDigits>( value, base ) );
        for (Size i=0; i<check.size(); ++i)
            CIE_TEST_CHECK( test[i] == check[i] );
    }

    {
        CIE_TEST_CASE_INIT( "base10" )

        Array value     = { 1, 0, 1, 1 };
        base            = 2;
        NT check        = 13;
        NT test;

        CIE_TEST_REQUIRE_NOTHROW( test = base10( value, base ) );
        CIE_TEST_CHECK( test == check );

        base            = 5;
        check           = 151;
        CIE_TEST_REQUIRE_NOTHROW( test = base10( value, base ) );
        CIE_TEST_CHECK( test == check );
    }

    {
        CIE_TEST_CASE_INIT( "reinterpretBase" )

        Array value     = { 1, 0, 1, 1 };
        NT base         = 2;
        NT check        = 13;
        NT test;

        CIE_TEST_REQUIRE_NOTHROW( test = reinterpretBase( value, base ) );
        CIE_TEST_CHECK( test == check );

        base            = 5;
        check           = 151;
        CIE_TEST_REQUIRE_NOTHROW( test = reinterpretBase( value, base ) );
        CIE_TEST_CHECK( test == check );
    }

} // CIE_TEST_CASE baseN


} // namespace cie::utils