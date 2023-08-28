// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/types/inc/vectortypes.hpp"
#include "packages/types/inc/typeoperations.hpp"

namespace cie::linalg {

CIE_TEST_CASE( "Norms and normalization", "[utilities]" )
{
    DoubleVector vector = {0.0, 0.0, 0.0, 0.0, 0.0};
    CIE_TEST_CHECK( norm<DoubleVector>(vector) == Approx(0.0) );
    CIE_TEST_REQUIRE_THROWS( normalize(vector) );

    vector = { -1.0, -4.0, 8.0 };
    CIE_TEST_CHECK( norm<DoubleVector>(vector) == Approx(9.0) );
    CIE_TEST_REQUIRE_NOTHROW( normalize(vector) );
    CIE_TEST_CHECK( vector[0] == Approx( -1/9.0 ) );
    CIE_TEST_CHECK( vector[1] == Approx( -4/9.0 ) );
    CIE_TEST_CHECK( vector[2] == Approx( 8/9.0 ) );
}

} // namespace cie::linalg