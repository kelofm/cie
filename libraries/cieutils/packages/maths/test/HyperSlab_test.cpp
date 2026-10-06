// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/maths/inc/HyperSlab.hpp"

// --- STL Includes ---
#include <algorithm>
#include <array>


namespace cie {


template <unsigned D, class T>
void TestHyperSlab() {
    {
        std::array<T,2> data{2, 5};
        HyperSlab<1,T> slab(data);
        std::array<T,1> range{8};
        std::array<T,3> output{};

        CIE_TEST_CHECK(slab.expand(range, output));
        CIE_TEST_CHECK(output[0] == 2);
        CIE_TEST_CHECK(output[1] == 3);
        CIE_TEST_CHECK(output[2] == 4);
    }

    if constexpr (D == 2) {
        std::array<T,4> data{1, 3, 2, 4};
        HyperSlab<2,T> slab(data);
        std::array<T,2> range{5, 6};
        std::array<T,4> output{};

        CIE_TEST_CHECK(slab.expand(range, output));
        CIE_TEST_CHECK(output[0] == 11);
        CIE_TEST_CHECK(output[1] == 12);
        CIE_TEST_CHECK(output[2] == 16);
        CIE_TEST_CHECK(output[3] == 17);
    }

    {
        std::array<T,2*D> data{};
        HyperSlab<D,T> slab(data);
        std::array<T,D> range{};
        range.fill(4);
        std::array<T,0> output{};

        CIE_TEST_CHECK(slab.empty());
        CIE_TEST_CHECK(slab.count() == 0);
        CIE_TEST_CHECK(slab.expand(range, output));
    }
} // void TestHyperSlab


CIE_TEST_CASE("Tiling", "[maths]") {
    CIE_TEST_CASE_INIT("Tiling")

    using T = unsigned;
    const HyperSlab<2,T> target(std::array<T,4>{1, 6, 2, 7});
    const Tiling<2,T> tiling(std::array<T,2>{2, 3});
    std::array<HyperSlab<2,T>,6> output;

    CIE_TEST_CHECK(tiling.count(target) == output.size());
    CIE_TEST_CHECK(tiling.tile(target, output));
    CIE_TEST_CHECK(std::equal(output[0].data().begin(), output[0].data().end(), std::array<T,4>{1, 3, 2, 5}.begin()));
    CIE_TEST_CHECK(std::equal(output[1].data().begin(), output[1].data().end(), std::array<T,4>{3, 5, 2, 5}.begin()));
    CIE_TEST_CHECK(std::equal(output[2].data().begin(), output[2].data().end(), std::array<T,4>{5, 6, 2, 5}.begin()));
    CIE_TEST_CHECK(std::equal(output[3].data().begin(), output[3].data().end(), std::array<T,4>{1, 3, 5, 7}.begin()));
    CIE_TEST_CHECK(std::equal(output[4].data().begin(), output[4].data().end(), std::array<T,4>{3, 5, 5, 7}.begin()));
    CIE_TEST_CHECK(std::equal(output[5].data().begin(), output[5].data().end(), std::array<T,4>{5, 6, 5, 7}.begin()));

    std::array<HyperSlab<2,T>,5> undersizedOutput;
    CIE_TEST_CHECK(!tiling.tile(target, undersizedOutput));
} // CIE_TEST_CASE Tiling


CIE_TEST_CASE("HyperSlab", "[maths]") {
    CIE_TEST_CASE_INIT("HyperSlab")

    {
        CIE_TEST_CASE_INIT("unsigned short")
        using T = unsigned short;

        {
            CIE_TEST_CASE_INIT("1D")
            TestHyperSlab<1,T>();
        }

        {
            CIE_TEST_CASE_INIT("2D")
            TestHyperSlab<2,T>();
        }

        {
            CIE_TEST_CASE_INIT("3D")
            TestHyperSlab<3,T>();
        }
    }

    {
        CIE_TEST_CASE_INIT("std::size_t")
        using T = std::size_t;

        {
            CIE_TEST_CASE_INIT("1D")
            TestHyperSlab<1,T>();
        }

        {
            CIE_TEST_CASE_INIT("2D")
            TestHyperSlab<2,T>();
        }

        {
            CIE_TEST_CASE_INIT("3D")
            TestHyperSlab<3,T>();
        }
    }
} // CIE_TEST_CASE HyperSlab


} // namespace cie
