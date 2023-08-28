// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Unique.hpp"


namespace cie::ct {


CIE_TEST_CASE("IsUnique", "[compile_time]")
{
    CIE_TEST_CASE_INIT("IsUnique")

    CIE_TEST_CHECK(IsUnique<>);
    CIE_TEST_CHECK(IsUnique<int>);
    CIE_TEST_CHECK(IsUnique<int,float>);
    CIE_TEST_CHECK(!IsUnique<int,int>);
    CIE_TEST_CHECK(!IsUnique<float,int,int>);
    CIE_TEST_CHECK(!IsUnique<int,float,int>);
    CIE_TEST_CHECK(!IsUnique<int,int,float>);
}


} // namespace cie::ct