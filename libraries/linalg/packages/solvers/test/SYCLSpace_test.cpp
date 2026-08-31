#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/solvers/inc/SYCLCSROperator.hpp"
#include "packages/solvers/inc/SYCLMaskedCSROperator.hpp"

// --- STL Includes ---
#include <vector>
#include <numeric>


namespace cie::linalg {


CIE_TEST_CASE("SYCLSpace", "[solvers]") {
    CIE_TEST_CASE_INIT("SYCLSpace")
    using T = double;

    std::shared_ptr<SYCLSpace<T>> pSpace;
    auto pQueue = std::make_shared<sycl::queue>(sycl::default_selector_v);
    CIE_TEST_REQUIRE_NOTHROW(pSpace = std::make_shared<SYCLSpace<T>>(pQueue));

    {
        CIE_TEST_CASE_INIT("innerProduct")

        for (const std::size_t size : {1e1+1, 1e2+2, 1e3+3, 1e4+4, 1e5+5, 1e6+6, 1e7+7}) {
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

            T product = 0;
            CIE_TEST_REQUIRE_NOTHROW(product = pSpace->innerProduct(left, right));
            CIE_TEST_CHECK(product == Approx(-T(size) * T(size - 1) / T(2)));
        }
    } // innerProduct

    {
        CIE_TEST_CASE_INIT("assign (from SYCL to SYCL)")
        for (std::size_t size : {1e1, 1e2, 1e3, 1e4, 1e5}) {
            std::vector<T> source(size);
            std::iota(source.begin(), source.end(), T(0));
            std::vector<T> target(source.size());

            auto deviceSource = pSpace->makeVector(source.size());
            auto deviceTarget = pSpace->makeVector(target.size());

            pSpace->assign(deviceSource, source);
            pSpace->assign(deviceTarget, deviceSource);
            pSpace->assign(target, deviceTarget);

            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent) {
                CIE_TEST_CHECK(target[iComponent] == source[iComponent]);
                CIE_TEST_CHECK(target[iComponent] == T(iComponent));
            }
        }
    }

    {
        CIE_TEST_CASE_INIT("add")
        for (std::size_t size : {1e1, 1e2, 1e3, 1e4, 1e5}) {
            std::vector<T> source(size);
            std::iota(source.begin(), source.end(), T(0));
            std::vector<T> target(source.size());
            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent)
                target[iComponent] = 1e1 * source[iComponent];

            auto deviceSource = pSpace->makeVector(source.size());
            auto deviceTarget = pSpace->makeVector(target.size());

            pSpace->assign(deviceSource, source);
            pSpace->assign(deviceTarget, target);
            pSpace->add(deviceTarget, deviceSource, 1e2);
            pSpace->assign(target, deviceTarget);

            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent) {
                CIE_TEST_CHECK(target[iComponent] == Approx((1e1 + 1e2) * source[iComponent]));
            }

            pSpace->add(deviceTarget, deviceSource, 1);
            pSpace->assign(target, deviceTarget);

            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent) {
                CIE_TEST_CHECK(target[iComponent] == Approx((1e1 + 1e2 + 1) * source[iComponent]));
            }
        }
    }

    {
        CIE_TEST_CASE_INIT("scale (uniform)")
        for (std::size_t size : {1e1, 1e2, 1e3, 1e4, 1e5}) {
            std::vector<T> target(size);
            std::iota(target.begin(), target.end(), T(0));
            auto deviceTarget = pSpace->makeVector(target.size());

            pSpace->assign(deviceTarget, target);
            pSpace->scale(deviceTarget, 2);
            pSpace->assign(target, deviceTarget);

            for (std::size_t iComponent=0ul; iComponent<target.size(); ++iComponent) {
                CIE_TEST_CHECK(target[iComponent] == Approx(T(iComponent * 2)));
            }
        }
    }

    {
        CIE_TEST_CASE_INIT("scale (component-wise)")
        for (std::size_t size : {1e1, 1e2, 1e3, 1e4, 1e5}) {
            std::vector<T> source(size);
            std::iota(source.begin(), source.end(), T(0));
            std::vector<T> target(source.size());
            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent)
                target[iComponent] = 1e1 * source[iComponent];

            auto deviceSource = pSpace->makeVector(source.size());
            auto deviceTarget = pSpace->makeVector(target.size());

            pSpace->assign(deviceSource, source);
            pSpace->assign(deviceTarget, target);
            pSpace->scale(deviceTarget, deviceSource, 1e2);
            pSpace->assign(target, deviceTarget);

            for (std::size_t iComponent=0ul; iComponent<source.size(); ++iComponent) {
                CIE_TEST_CHECK(target[iComponent] == Approx(1e1 * 1e2 * source[iComponent] * source[iComponent]));
            }
        }
    }

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

    {
        CIE_TEST_CASE_INIT("SpMV (masked)")
        // Define a matrix.
        const std::vector<T> entries {1, 2, 100, 3, 4};
        const std::vector<int> columnIndices {0, 2, 1, 0, 2};
        const std::vector<int> rowExtents {0, 2, 3, 5};
        const std::vector<std::uint16_t> mask {1, 2, 1};
        const std::uint16_t threshold = 2;

        // Allocate memory for the matrix on the device.
        SYCLSpace<int> indexSpace(pQueue);
        auto deviceEntries = pSpace->makeVector(entries.size());
        auto deviceColumnIndices = indexSpace.makeVector(columnIndices.size());
        auto deviceRowExtents = indexSpace.makeVector(rowExtents.size());

        auto pMaskSpace = std::make_shared<SYCLSpace<std::uint16_t>>(pQueue);
        auto deviceMask = pMaskSpace->makeVector(3);
        pMaskSpace->assign(deviceMask, mask);

        // Copy the matrix from the host to the device.
        pSpace->assign(deviceEntries, entries);
        indexSpace.assign(deviceColumnIndices, columnIndices);
        indexSpace.assign(deviceRowExtents, rowExtents);

        // Construct a view over the device matrix.
        const CSRView<const T,const int> matrix(
            /*columnCount=*/3,
            {deviceRowExtents.get(), deviceRowExtents.size()},
            {deviceColumnIndices.get(), deviceColumnIndices.size()},
            {deviceEntries.get(), deviceEntries.size()});

        // Define the linear operator.
        SYCLMaskedCSROperator<int,T,std::uint16_t> op(
            matrix,
            deviceMask,
            threshold,
            pSpace,
            pMaskSpace);

        // Define the right hand and result vectors.
        const std::vector<T> rhs {5, 100, 6};
        std::vector<T> result {7, 1000, 8};
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
        CIE_TEST_REQUIRE(result.size() == 3);
        CIE_TEST_CHECK(result[0] == Approx(233.0));
        CIE_TEST_CHECK(result[1] == Approx(9 * 1000.0));
        CIE_TEST_CHECK(result[2] == Approx(462.0));
    } // SpMV
} // CIE_TEST_CASE "SYCLSpace"


} // namespace cie::linalg

#endif
