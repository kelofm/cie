#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLSpace.hpp"

// --- STL Includes ---
#include <vector>
#include <numeric>


namespace cie::linalg {


CIE_TEST_CASE("SYCLSpace", "[solvers]") {
    CIE_TEST_CASE_INIT("SYCLSpace")
    using T = float;

    std::shared_ptr<SYCLSpace<T>> pSpace;
    CIE_TEST_REQUIRE_NOTHROW(pSpace = std::make_shared<SYCLSpace<T>>(
        std::make_shared<sycl::queue>(sycl::default_selector_v)));

    const std::size_t size = 1e4;

    auto left = pSpace->makeVector(size);
    auto right = pSpace->makeVector(size);
    auto result = pSpace->makeVector(size);

    {
        std::vector<T> leftSource(size), rightSource(size);
        std::iota(leftSource.begin(), leftSource.end(), T(0));
        std::fill(rightSource.begin(), rightSource.end(), T(-1));
        pSpace->assign(left, leftSource);
        pSpace->assign(right, rightSource);
    }

    {
        CIE_TEST_CASE_INIT("innerProduct")
        T result = 0;
        CIE_TEST_REQUIRE_NOTHROW(result = pSpace->innerProduct(left, right));
        CIE_TEST_CHECK(result == Approx(-T(size) * (size - 1) / 2));
    }

} // CIE_TEST_CASE "SYCLSpace"


} // namespace cie::linalg

#endif
