// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/texture/inc/Image.hpp"
#include "packages/utility/inc/paths.hpp"


namespace cie::gl {


CIE_TEST_CASE( "Image", "[texture]" )
{
    CIE_TEST_CASE_INIT( "Image" )

    CIE_TEST_REQUIRE_NOTHROW( Image( GL_TEXTURE_SOURCE_PATH / "cage.png" ) );
    Image image( GL_TEXTURE_SOURCE_PATH / "cage.png" );

    CIE_TEST_CHECK( image.width() == 511 );
    CIE_TEST_CHECK( image.height() == 543 );
    CIE_TEST_CHECK( image.numberOfChannels() == 3 );
    CIE_TEST_CHECK( image.size() == 511 * 543 * 3 );

    CIE_TEST_CHECK_NOTHROW( image.resize(30, 30) );
    CIE_TEST_CHECK( image.width() == 30 );
    CIE_TEST_CHECK( image.height() == 30 );
    CIE_TEST_CHECK( image.numberOfChannels() == 3 );
    CIE_TEST_CHECK( image.size() == 30 * 30 * 3 );

    CIE_TEST_CHECK_NOTHROW( image.write(getOutputPath()/ "cage.png") );
    CIE_TEST_CHECK_NOTHROW( image.write(getOutputPath()/ "cage.jpg") );
    CIE_TEST_CHECK_NOTHROW( image.write(getOutputPath()/ "cage.tga") );
    CIE_TEST_CHECK_NOTHROW( image.write(getOutputPath()/ "cage.bmp") );

    CIE_TEST_CHECK_NOTHROW( image.clear() );
    CIE_TEST_CHECK( image.width() == 0 );
    CIE_TEST_CHECK( image.height() == 0 );
    CIE_TEST_CHECK( image.numberOfChannels() == 0 );
    CIE_TEST_CHECK( image.size() == 0 );

    CIE_TEST_CHECK_NOTHROW( image.load(getOutputPath()/ "cage.png") );
    CIE_TEST_CHECK_NOTHROW( image.load(getOutputPath()/ "cage.jpg") );
    CIE_TEST_CHECK_NOTHROW( image.load(getOutputPath()/ "cage.tga") );
    CIE_TEST_CHECK_NOTHROW( image.load(getOutputPath()/ "cage.bmp") );
    CIE_TEST_CHECK( image.width() == 30 );
    CIE_TEST_CHECK( image.height() == 30 );
    CIE_TEST_CHECK( image.numberOfChannels() == 3 );
    CIE_TEST_CHECK( image.size() == 30 * 30 * 3 );

    CIE_TEST_CHECK_NOTHROW( Image(image) );

    Image image2;
    CIE_TEST_CHECK_NOTHROW( image2 = image );
    CIE_TEST_CHECK( image2.width() == 30 );
    CIE_TEST_CHECK( image2.height() == 30 );
    CIE_TEST_CHECK( image2.numberOfChannels() == 3 );
    CIE_TEST_CHECK( image2.size() == 30 * 30 * 3 );
}


} // namespace cie::gl