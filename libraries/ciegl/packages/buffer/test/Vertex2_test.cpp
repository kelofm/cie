// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex2.hpp"


namespace cie::gl {


CIE_TEST_CASE( "Vertex2", "[buffer]" )
{
    CIE_TEST_CASE_INIT( "Vertex2" )

    auto p_attributes = AttributeContainer::SharedPointer( new AttributeContainer );

    CIE_TEST_REQUIRE_NOTHROW( makeVertex<Vertex2>( p_attributes ) );

    Vertex2::SharedPointer p_vertex;
    CIE_TEST_REQUIRE_NOTHROW( p_vertex = makeVertex<Vertex2>( p_attributes ) );

    CIE_TEST_CHECK_NOTHROW( p_vertex->x() = 1.0 );
    CIE_TEST_CHECK_NOTHROW( p_vertex->y() = 2.0 );
    CIE_TEST_CHECK( p_vertex->x() == Approx(1.0) );
    CIE_TEST_CHECK( p_vertex->y() == Approx(2.0) );

    CIE_TEST_REQUIRE_NOTHROW( p_vertex->position() );
    Size counter = 0;
    for ( const auto& component : p_vertex->position() )
        CIE_TEST_CHECK( component == Approx( ++counter ) );

    for ( auto& component : p_vertex->position() )
    {
        double tmp = (counter+1) / double(counter+2);
        component = tmp;
        counter += 2;
    }

    CIE_TEST_CHECK( p_vertex->x() == Approx(3.0 / 4.0) );
    CIE_TEST_CHECK( p_vertex->y() == Approx(5.0 / 6.0) );
}


} // namespace cie::gl