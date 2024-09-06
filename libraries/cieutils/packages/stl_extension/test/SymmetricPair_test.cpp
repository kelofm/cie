// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/SymmetricPair.hpp"

// --- STL Includes ---
#include <unordered_set>


namespace cie::utils {


CIE_TEST_CASE("SymmetricPair", "[stl_extension]")
{
    CIE_TEST_CASE_INIT("SymmetricPair")
    using P = SymmetricPair<int>;
    std::unordered_set<P,Hash<P>> set;

    P null {0, 0};
    set.insert(null);
    CIE_TEST_REQUIRE(set.size() == 1);

    CIE_TEST_CHECK(null == null);
    CIE_TEST_CHECK(!(null != null));
    CIE_TEST_CHECK(*set.begin() == null);

    P pair01 {0, 1};
    set.insert(pair01);
    CIE_TEST_CHECK(!(null == pair01));
    CIE_TEST_CHECK(null != pair01);
    CIE_TEST_CHECK(set.size() == 2);

    P pair10 {1, 0};
    set.insert(pair10);
    CIE_TEST_CHECK(!(null == pair10));
    CIE_TEST_CHECK(null != pair10);
    CIE_TEST_CHECK(pair10 == pair01);
    CIE_TEST_CHECK(!(pair10 != pair01));
    CIE_TEST_CHECK(set.size() == 2);
}


} // namespace cie::utils
