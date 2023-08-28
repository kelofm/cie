// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Select.hpp"


namespace cie::ct {



CIE_TEST_CASE("Select", "[compile_time]")
{

    CIE_TEST_CASE_INIT("Select")

    CIE_TEST_CHECK(std::is_same_v<Select<void,char,int>::Front, void>);
    CIE_TEST_CHECK(std::is_same_v<Select<void,char,int>::At<0>, void>);
    CIE_TEST_CHECK(std::is_same_v<Select<void,char,int>::At<1>, char>);
    CIE_TEST_CHECK(std::is_same_v<Select<void,char,int>::At<2>, int>);
    CIE_TEST_CHECK(std::is_same_v<Select<void,char,int>::Back, int>);

    CIE_TEST_CHECK(std::is_same_v<Select<void>::Front, void>);
    CIE_TEST_CHECK(std::is_same_v<Select<void>::At<0>, void>);
    CIE_TEST_CHECK(std::is_same_v<Select<void>::Back, void>);
}


} // namespace cie::ct