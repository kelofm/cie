// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/state_iterator.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <vector>
#include <array>
#include <iostream>


namespace cie::utils {


CIE_TEST_CASE( "StateIterator", "[stl_extension]" )
{
    typedef std::vector<int>                            SubContainer;
    typedef StaticArray<std::vector<int>,2>              Container;

    const Container container =
    {
        SubContainer({ 0, 1, 2 }),
        SubContainer({ 3, 4 })
    };

    CIE_TEST_REQUIRE_NOTHROW( makeStateIterator( container.begin(), container.end() ) );
    CIE_TEST_REQUIRE_NOTHROW( makeStateIterator( container ) );

    auto state = makeStateIterator(container);
    //StateIterator<typename Container::const_iterator> state(container);
    CIE_TEST_REQUIRE_NOTHROW( *state );
    CIE_TEST_REQUIRE_NOTHROW( ++state );
    CIE_TEST_REQUIRE_NOTHROW( state.reset() );

    // state 0 -> 0,3
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 0 );     CIE_TEST_CHECK( *(*state)[1] == 3 );
    CIE_TEST_CHECK_NOTHROW( ++state );

    // state 0 -> 1,3
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 1 );     CIE_TEST_CHECK( *(*state)[1] == 3 );
    CIE_TEST_CHECK_NOTHROW( ++state );

    // state 0 -> 2,3
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 2 );     CIE_TEST_CHECK( *(*state)[1] == 3 );
    CIE_TEST_CHECK_NOTHROW( ++state );

    // state 0 -> 0,4
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 0 );     CIE_TEST_CHECK( *(*state)[1] == 4 );
    CIE_TEST_CHECK_NOTHROW( ++state );

    // state 0 -> 1,4
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 1 );     CIE_TEST_CHECK( *(*state)[1] == 4 );
    CIE_TEST_CHECK_NOTHROW( ++state );

    // state 0 -> 2,4
    CIE_TEST_REQUIRE( (*state).size() == container.size() );
    CIE_TEST_CHECK( *(*state)[0] == 2 );     CIE_TEST_CHECK( *(*state)[1] == 4 );
    CIE_TEST_CHECK_NOTHROW( ++state );

} // CIE_TEST_CASE( "StateIterator")


CIE_TEST_CASE( "InternalStateIterator", "[stl_extension]" )
{
    const Size dimension        = 2;
    StaticArray<int,3> container = { 0, 1, 2 };
    auto state                  = makeInternalStateIterator( container, dimension );

    for (Size i=0; i<container.size(); ++i)
        for (Size j=0; j<container.size(); ++j)
        {
            CIE_TEST_REQUIRE( (*state).size() == dimension );
            CIE_TEST_CHECK( *(*state)[0] == container[j] );     CIE_TEST_CHECK( *(*state)[1] == container[i] );
            CIE_TEST_CHECK_NOTHROW( ++state );
        }
} // CIE_TEST_CASE( "InternalStateIterator")


} // namespace cie::utils