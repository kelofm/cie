// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/file/inc/STLPart.hpp"
#include "packages/utility/inc/paths.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::gl {


CIE_TEST_CASE( "STLPart", "[file]" )
{
    CIE_TEST_CASE_INIT( "STLPart" )

    PartPtr p_part;

    CIE_TEST_REQUIRE_NOTHROW( p_part.reset(
        new STLPart( GL_GEOMETRY_SOURCE_PATH / "stl/Bunny.stl" ))
    );
}


} // namespace cie::gl