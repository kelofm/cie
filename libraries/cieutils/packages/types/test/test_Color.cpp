// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/types/inc/Color.hpp"


namespace cie {


CIE_TEST_CASE("Color", "[utility]")
{
    CIE_TEST_CASE_INIT("Color")

    CIE_TEST_CHECK(RGBAColor::Black.hex() == 0x0);
    CIE_TEST_CHECK(RGBAColor::White.hex() == 0xffffff);
    CIE_TEST_CHECK(RGBAColor::TUMBlue.hex() == 0x5293);
}


} // namespace cie