// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/AttributeContainer.hpp"
#include "packages/buffer/inc/AbsVertex.hpp"


namespace cie::gl {


CIE_TEST_CASE( "AttributeContainer", "[buffer]" )
{
    CIE_TEST_CASE_INIT( "AttributeContainer" )

    CIE_TEST_REQUIRE_NOTHROW(AttributeContainer::SharedPointer( new AttributeContainer ) );

    AttributeContainer::SharedPointer attributeContainer(
        new AttributeContainer
    );

    // TODO
}


} // namespace cie::gl