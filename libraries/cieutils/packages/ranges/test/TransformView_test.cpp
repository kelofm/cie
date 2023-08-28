// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/ranges/inc/TransformView.hpp"

// --- STL Includes ---
#include <deque>


namespace cie::utils {


CIE_TEST_CASE( "TransformView", "[ranges]" )
{
    CIE_TEST_CASE_INIT( "TransformView" )

    using ValueType     = float;
    using ContainerType = std::deque<ValueType>;

    ContainerType container { 0, 1, 2, 3, 4 };

    {
        using OutputType    = int;

        CIE_TEST_CHECK_NOTHROW(
            makeTransformView<OutputType>( container,
                                           [](auto value) -> OutputType {return value;} )
        );

        auto container_view = makeTransformView<OutputType>( container,
                                                             [](auto value) -> OutputType {return value;} );

        for ( Size index=0; index<container.size(); ++index )
            CIE_TEST_CHECK( container[index] == container_view[index] );
    }
}


} // namespace cie::utils