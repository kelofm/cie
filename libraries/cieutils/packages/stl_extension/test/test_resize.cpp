// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/resize.hpp"

// --- STL Includes ---
#include <vector>
#include <array>


namespace cie {


CIE_TEST_CASE( "Concept utilities", "[concepts]" )
{
    std::vector<int>    vector(1);
    std::array<int,3>   array;

    CIE_TEST_CHECK_NOTHROW( utils::resize(vector,3) );
    CIE_TEST_CHECK( vector.size() == 3 );

    CIE_TEST_CHECK_NOTHROW( utils::resize(array, 3) );
    CIE_TEST_CHECK( array.size() == 3 );
}


}
