#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/solvers/inc/SYCLCSROperator.hpp"

// --- STL Includes ---
#include <vector>
#include <numeric>


namespace cie::linalg {


CIE_TEST_CASE("SYCLSpace", "[solvers]") {
    CIE_TEST_CASE_INIT("SYCLSpace")
    using T = float;

    std::shared_ptr<SYCLSpace<T>> pSpace;
    auto pQueue = std::make_shared<sycl::queue>(sycl::default_selector_v);
    CIE_TEST_REQUIRE_NOTHROW(pSpace = std::make_shared<SYCLSpace<T>>(pQueue));

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
    } // innerProduct

    {
        CIE_TEST_CASE_INIT("SpMV")
        // Define a matrix.
        const std::vector<T> entries {1, 2, 3, 4};
        const std::vector<int> columnIndices {0, 1, 0, 1};
        const std::vector<int> rowExtents {0, 2, 4};

        // Allocate memory for the matrix on the device.
        SYCLSpace<int> indexSpace(pQueue);
        auto deviceEntries = pSpace->makeVector(entries.size());
        auto deviceColumnIndices = indexSpace.makeVector(columnIndices.size());
        auto deviceRowExtents = indexSpace.makeVector(rowExtents.size());

        // Copy the matrix from the host to the device.
        pSpace->assign(deviceEntries, entries);
        indexSpace.assign(deviceColumnIndices, columnIndices);
        indexSpace.assign(deviceRowExtents, rowExtents);

        // Construct a view over the device matrix.
        const CSRView<const T,const int> matrix(
            /*columnCount=*/2,
            {deviceRowExtents.get(), deviceRowExtents.size()},
            {deviceColumnIndices.get(), deviceColumnIndices.size()},
            {deviceEntries.get(), deviceEntries.size()});

        // Define the linear operator.
        SYCLCSROperator<int,T> op(matrix, pSpace);

        // Define the right hand and result vectors.
        const std::vector<T> rhs {5, 6};
        std::vector<T> result {7, 8};
        auto deviceRHS = pSpace->makeVector(matrix.rowCount());
        auto deviceResult = pSpace->makeVector(matrix.columnCount());
        pSpace->assign(deviceRHS, rhs);
        pSpace->assign(deviceResult, result);

        // Perform the SpMV.
        CIE_TEST_REQUIRE_NOTHROW(op.product(
            9,
            deviceRHS,
            10,
            deviceResult));
        pSpace->assign(result, deviceResult);
        CIE_TEST_REQUIRE(result.size() == 2);
        CIE_TEST_CHECK(result[0] == Approx(233.0));
        CIE_TEST_CHECK(result[1] == Approx(462.0));
    } // SpMV
} // CIE_TEST_CASE "SYCLSpace"


} // namespace cie::linalg

#endif
