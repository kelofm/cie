// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/OuterProduct.hpp"

// --- STL Includes ---
#include <array>
#include <span>


namespace cie::maths {


CIE_TEST_CASE("OuterProduct", "[maths]")
{
    CIE_TEST_CASE_INIT("next in HyperSlab")

    const std::array<unsigned,4> slabData{1, 4, 2, 5};
    const HyperSlab<2,unsigned> slab(slabData);
    std::array<unsigned,2> state{2, 3};
    std::span<unsigned,2> stateSpan(state);

    CIE_TEST_CHECK(OuterProduct<2>::next(slab, stateSpan));
    CIE_TEST_CHECK(state == std::array<unsigned,2>{3, 3});

    CIE_TEST_CHECK(OuterProduct<2>::next(slab, stateSpan));
    CIE_TEST_CHECK(state == std::array<unsigned,2>{1, 4});

    CIE_TEST_CHECK(OuterProduct<2>::next(slab, stateSpan));
    CIE_TEST_CHECK(state == std::array<unsigned,2>{2, 4});

    state = {3, 4};
    CIE_TEST_CHECK(!OuterProduct<2>::next(slab, stateSpan));
    CIE_TEST_CHECK(state == std::array<unsigned,2>{1, 2});
}


} // namespace cie::maths
