// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/MarchingContainer.hpp"


namespace cie::utils {


CIE_TEST_CASE( "MarchingContainer", "[stl_extension]" )
{
    CIE_TEST_CASE_INIT( "MarchingContainer" )

    using ValueType     = int;
    using ContainerType = MarchingContainer<ValueType>;

    for ( Size capacity=3; capacity<7; ++capacity )
    {
        CIE_TEST_CASE_INIT( "capacity " + std::to_string(capacity) )

        {
            CIE_TEST_REQUIRE_NOTHROW( ContainerType(capacity) );
            ContainerType container( capacity );

            ValueType counter = 0;
            CIE_TEST_CHECK_NOTHROW( container.push_back( ++counter ) );
            CIE_TEST_CHECK_NOTHROW( container.emplace_back( ++counter ) );

            CIE_TEST_CHECK_NOTHROW( container.back() );
            CIE_TEST_CHECK( container.back() == counter );

            CIE_TEST_CHECK_NOTHROW( container.back() = ++counter );
            CIE_TEST_CHECK( container.back() == counter );
            auto tmp = container.back();

            CIE_TEST_CHECK_NOTHROW( container[0] );
            CIE_TEST_CHECK( container[0] == 1 );

            CIE_TEST_CHECK_NOTHROW( container[0] = ++counter );
            CIE_TEST_CHECK( container[0] == counter );

            CIE_TEST_CHECK_NOTHROW( container[1] );
            CIE_TEST_CHECK( container[1] == tmp );

            CIE_TEST_CHECK_NOTHROW( container[1] = ++counter );
            CIE_TEST_CHECK( container[1] == counter );

            CIE_TEST_CHECK_NOTHROW( container.size() );
            CIE_TEST_CHECK( container.size() == 2 );

            CIE_TEST_CHECK_NOTHROW( container[ container.size() - 1 ] );
            CIE_TEST_CHECK( container[ container.size()-1 ] == container.back() );

            #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
            CIE_TEST_CHECK_THROWS( container[container.size()] );
            #endif

            for ( Size i=0; i<3*capacity; ++i )
                CIE_TEST_CHECK_NOTHROW( container.push_back( ValueType(++counter) ) );

            #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
            CIE_TEST_CHECK_THROWS( container[ container.size() - capacity - 1 ] );
            #endif

            for ( Size i=0; i<capacity; ++i )
            {
                Size index = container.size() - capacity + i;
                CIE_TEST_REQUIRE_NOTHROW( container[index] );
                CIE_TEST_CHECK( container[index] == ValueType(counter-capacity+i+1) );
            }

            CIE_TEST_CHECK( container.back() == counter );

            #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
            CIE_TEST_CHECK_THROWS( container[container.size()] );
            #endif
        }

        {
            ValueType initializer = 12;
            ContainerType container( capacity, initializer );

            CIE_TEST_CHECK( container.size() == capacity );

            for ( Size i=0; i<capacity; ++i )
            {
                CIE_TEST_CHECK_NOTHROW( container[i] );
                CIE_TEST_CHECK( container[i] == initializer );
            }

            for ( Size i=capacity; i<4*capacity+1; ++i )
                CIE_TEST_CHECK_NOTHROW( container.emplace_back( ValueType(i) ) );

            for ( Size i=container.size()-capacity; i<container.size(); ++i )
            {
                CIE_TEST_CHECK_NOTHROW( container[i] );
                CIE_TEST_CHECK( container[i] == ValueType(i) );
            }
        }
    }
}


} // namespace cie::utils