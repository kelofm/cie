// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex3.hpp"


namespace cie::gl {


CIE_TEST_CASE( "Vertex3", "[buffer]" )
{
    CIE_TEST_CASE_INIT( "Vertex3" )

    auto p_attributes = AttributeContainer::SharedPointer( new AttributeContainer );

    CIE_TEST_REQUIRE_NOTHROW( makeVertex<Vertex3>( p_attributes ) );

    Vertex3::SharedPointer p_vertex;
    CIE_TEST_REQUIRE_NOTHROW( p_vertex = makeVertex<Vertex3>( p_attributes ) );

    CIE_TEST_CHECK_NOTHROW( p_vertex->x() = 1.0 );
    CIE_TEST_CHECK_NOTHROW( p_vertex->y() = 2.0 );
    CIE_TEST_CHECK_NOTHROW( p_vertex->z() = 3.0 );
    CIE_TEST_CHECK( p_vertex->x() == Approx(1.0) );
    CIE_TEST_CHECK( p_vertex->y() == Approx(2.0) );
    CIE_TEST_CHECK( p_vertex->z() == Approx(3.0) );

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

    CIE_TEST_CHECK( p_vertex->x() == Approx(4.0 / 5.0) );
    CIE_TEST_CHECK( p_vertex->y() == Approx(6.0 / 7.0) );
    CIE_TEST_CHECK( p_vertex->z() == Approx(8.0 / 9.0) );
}


} // namespace cie::gl