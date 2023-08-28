// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/types/inc/types.hpp"


namespace cie {


CIE_TEST_CASE( "intPow", "[maths]" )
{
    typedef Size    NT;

    NT base     = 2;
    NT check    = 1;
    for (Size i=0; i<17; ++i)
    {
        CIE_TEST_REQUIRE_NOTHROW( intPow(base,i) );
        CIE_TEST_CHECK( intPow(base,i) == check );
        check *= base;
    }

    base        = 13;
    check       = 1;
    for (Size i=0; i<5; ++i)
    {
        CIE_TEST_REQUIRE_NOTHROW( intPow(base,i) );
        CIE_TEST_CHECK( intPow(base,i) == check );
        check *= base;
    }
}


} // namespace cie