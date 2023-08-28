// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/concepts/inc/logical.hpp"


namespace cie {


template <concepts::Any<int,char,double> T>
bool anyTest(T t) {return true;}


template <concepts::None<int,char,double> T>
bool anyTest(T t) {return false;}


CIE_TEST_CASE("Any - None", "[concepts]")
{
    CIE_TEST_CASE_INIT("Any - None")

    const int i = 0;
    const char c = 0;
    const double d = 0;
    const float f = 0;

    CIE_TEST_CHECK(anyTest(i));
    CIE_TEST_CHECK(anyTest(c));
    CIE_TEST_CHECK(anyTest(d));
    CIE_TEST_CHECK(!anyTest(f));
}


} // namespace cie
