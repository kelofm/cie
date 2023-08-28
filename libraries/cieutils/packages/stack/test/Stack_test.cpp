// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stack/inc/Stack.hpp"


namespace cie::utils {


CIE_TEST_CASE( "Stack", "[stack]" )
{
    CIE_TEST_CASE_INIT( "Stack" )

    Stack<Size> stack( 0, 5 );

    for ( Size i=0; i<3; ++i )
        stack.push_back( i );

    CIE_TEST_CHECK( stack.back() == 2 );

    for ( Size i=0; i<3; ++i )
        CIE_TEST_CHECK( stack[i] == i );

    for ( Size i=3; i<12; ++i )
        stack.push_back(i);

    for ( Size i=7; i<12; ++i )
        CIE_TEST_CHECK( i == stack[i] );

    CIE_TEST_CHECK( stack.back() == 11 );
}


} // namespace cie::utils