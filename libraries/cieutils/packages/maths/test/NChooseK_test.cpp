// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/NChooseK.hpp"


namespace cie::utils {


CIE_TEST_CASE( "NChooseK", "[maths]" )
{
    CIE_TEST_CASE_INIT( "NChooseK" )

    // 4 choose 0
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(4, 0) );
        NChooseK test( 4, 0 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 0 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 1 );

        CIE_TEST_CHECK( !test.next() );
    }

    // 4 choose 1
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(4, 1) );
        NChooseK test( 4, 1 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 1 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 4 );
        
        CIE_TEST_CHECK( *test.begin() == 0 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *test.begin() == 1 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *test.begin() == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *test.begin() == 3 );

        CIE_TEST_CHECK( !test.next() );
    }

    // 4 choose 2
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(4, 2) );
        NChooseK test( 4, 2 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 2 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 6 );
        
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 1 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 2 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( !test.next() );
    }

    // 4 choose 3
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(4, 3) );
        NChooseK test( 4, 3 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 3 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 4 );
        
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 1 );
        CIE_TEST_CHECK( *(test.begin()+2) == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 1 );
        CIE_TEST_CHECK( *(test.begin()+2) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );
        CIE_TEST_CHECK( *(test.begin()+2) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );
        CIE_TEST_CHECK( *(test.begin()+2) == 3 );

        CIE_TEST_CHECK( !test.next() );
    }

    // 4 choose 4
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(4, 4) );
        NChooseK test( 4, 4 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 4 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 1 );
        
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 1 );
        CIE_TEST_CHECK( *(test.begin()+2) == 2 );
        CIE_TEST_CHECK( *(test.begin()+3) == 3 );

        CIE_TEST_CHECK( !test.next() );
    }

    // 5 choose 2
    {
        CIE_TEST_REQUIRE_NOTHROW( NChooseK(5, 2) );
        NChooseK test( 5, 2 );

        CIE_TEST_CHECK( std::distance(test.begin(), test.end()) == 2 );
        CIE_TEST_CHECK( test.numberOfPermutations() == 10 );
        
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 1 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 0 );
        CIE_TEST_CHECK( *(test.begin()+1) == 4 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 2 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 1 );
        CIE_TEST_CHECK( *(test.begin()+1) == 4 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 2 );
        CIE_TEST_CHECK( *(test.begin()+1) == 3 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 2 );
        CIE_TEST_CHECK( *(test.begin()+1) == 4 );

        CIE_TEST_CHECK( test.next() );
        CIE_TEST_CHECK( *(test.begin()) == 3 );
        CIE_TEST_CHECK( *(test.begin()+1) == 4 );

        CIE_TEST_CHECK( !test.next() );
    }
}


} // namespace cie::maths