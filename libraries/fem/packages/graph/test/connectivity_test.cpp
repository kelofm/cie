// --- Utility Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"


namespace cie::fem {


CIE_TEST_CASE("scanConnectivities", "[graph]")
{
    CIE_TEST_CASE_INIT("scanConnectivities")
    using Basis = maths::Polynomial<double>;

    DynamicArray<Basis> basisFunctions {
        Basis(Basis::Coefficients {0.5, -0.5}),
        Basis(Basis::Coefficients {0.5, 0.5})
    };

    StaticArray<double,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    DynamicArray<std::pair<BoundaryID,unsigned>> connectivities;
    const auto functor = [&connectivities] (BoundaryID boundaryID, unsigned iAnsatz) {
        connectivities.emplace_back(boundaryID, iAnsatz);
    };

    { // Linear shape functions on a square
        CIE_TEST_CASE_INIT("2D")
        maths::AnsatzSpace<Basis,2> ansatzSpace(basisFunctions);
        connectivities.clear();

        CIE_TEST_CHECK_NOTHROW(scanConnectivities(
            ansatzSpace,
            functor,
            samples.begin(),
            samples.end(),
            1e-10
        ));

        CIE_TEST_CHECK(connectivities.size() == 8);

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 0u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 2u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 1u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 3u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 0u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 1u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 2u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 3u)
        ) != connectivities.end());
    } // 2D

    { // Linear shape functions on a cube
        CIE_TEST_CASE_INIT("3D")
        maths::AnsatzSpace<Basis,3> ansatzSpace(basisFunctions);
        connectivities.clear();

        CIE_TEST_CHECK_NOTHROW(scanConnectivities(
            ansatzSpace,
            functor,
            samples.begin(),
            samples.end(),
            1e-10
        ));

        CIE_TEST_CHECK(connectivities.size() == 24);

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 0u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 0u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 0), 0u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 1u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 1u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 0), 1u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 2u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 2u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 0), 2u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 3u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 3u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 0), 3u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 4u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 4u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 1), 4u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 5u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 0), 5u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 1), 5u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 0), 6u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 6u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 1), 6u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(0, 1), 7u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(1, 1), 7u)
        ) != connectivities.end());

        CIE_TEST_CHECK(std::find(
            connectivities.begin(),
            connectivities.end(),
            std::make_pair(BoundaryID(2, 1), 7u)
        ) != connectivities.end());
    }
}



CIE_TEST_CASE("AnsatzMap", "[graph]")
{
    CIE_TEST_CASE_INIT("AnsatzMap")
    using Basis = maths::Polynomial<double>;

    const DynamicArray<Basis> basisFunctions {
        Basis(Basis::Coefficients { 0.5, -0.5}),
        Basis(Basis::Coefficients { 0.5,  0.5}),
        Basis(Basis::Coefficients { 1.0,  0.0, -1.0})
    };

    const StaticArray<double,5> samples {-1.0,
                                         -0.5,
                                          0.0,
                                          0.5,
                                          1.0};

    {
        CIE_TEST_CASE_INIT("1D")
        maths::AnsatzSpace<Basis,1> ansatzSpace(basisFunctions);
        const auto map = makeAnsatzMap(ansatzSpace,
                                       {samples.data(), samples.size()},
                                       utils::Comparison<double>(1e-10, 1e-10));

        {
            DynamicArray<std::pair<Size,Size>> coincidentPairs;
            const OrientedBoundary<1> first("+x", "+x");
            const OrientedBoundary<1> second("-x", "+x");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 1);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 1);
            CIE_TEST_CHECK(coincidentPairs.front().first == 1);
            CIE_TEST_CHECK(coincidentPairs.front().second == 1);
            coincidentPairs.clear();
        }

        {
            DynamicArray<std::pair<Size,Size>> coincidentPairs;
            const OrientedBoundary<1> first("+x", "-x");
            const OrientedBoundary<1> second("-x", "-x");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 1);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 1);
            CIE_TEST_CHECK(coincidentPairs.front().first == 0);
            CIE_TEST_CHECK(coincidentPairs.front().second == 0);
            coincidentPairs.clear();
        }
    } // 1D

    {
        CIE_TEST_CASE_INIT("2D")
        maths::AnsatzSpace<Basis,2> ansatzSpace(basisFunctions);
        const auto map = makeAnsatzMap(ansatzSpace,
                                       {samples.data(), samples.size()},
                                       utils::Comparison<double>(1e-10, 1e-10));
        DynamicArray<std::pair<Size,Size>> coincidentPairs;

        {
            OrientedBoundary<2> first("+x+y", "+x");
            OrientedBoundary<2> second("-x+y", "+x");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 3);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(1, 1)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(4, 4)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(7, 7)
            ) != coincidentPairs.end());
            coincidentPairs.clear();
        }

        {
            OrientedBoundary<2> first("-x+y", "+x");
            OrientedBoundary<2> second("+x+y", "+x");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 3);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(1, 1)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(4, 4)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(7, 7)
            ) != coincidentPairs.end());
            coincidentPairs.clear();
        }

        {
            OrientedBoundary<2> first("+x+y", "+y");
            OrientedBoundary<2> second("+x-y", "+y");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 3);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(3, 3)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(4, 4)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(5, 5)
            ) != coincidentPairs.end());
            coincidentPairs.clear();
        }

        {
            OrientedBoundary<2> first("+x-y", "+y");
            OrientedBoundary<2> second("+x+y", "+y");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 3);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(3, 3)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(4, 4)
            ) != coincidentPairs.end());

            CIE_TEST_CHECK(std::find(
                coincidentPairs.begin(),
                coincidentPairs.end(),
                std::pair<Size,Size>(5, 5)
            ) != coincidentPairs.end());
            coincidentPairs.clear();
        }
    } // 2D

    {
        CIE_TEST_CASE_INIT("3D")
        maths::AnsatzSpace<Basis,3> ansatzSpace(basisFunctions);
        const auto map = makeAnsatzMap(ansatzSpace,
                                       {samples.data(), samples.size()},
                                       utils::Comparison<double>(1e-10, 1e-8));
        DynamicArray<std::pair<Size,Size>> coincidentPairs;

        {
            OrientedBoundary<3> first("+x+y+z", "+x");
            OrientedBoundary<3> second("-x+y+z", "+x");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 9);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
            for (auto pair : StaticArray<std::pair<Size,Size>,9> {
                    { 1,  1},
                    { 4,  4},
                    {10, 10},
                    {13, 13},
                    { 7,  7},
                    {16, 16},
                    {19, 19},
                    {22, 22},
                    {25, 25}
                }) {
                CIE_TEST_CHECK(std::find(
                    coincidentPairs.begin(),
                    coincidentPairs.end(),
                    pair
                ) != coincidentPairs.end());
            }
            coincidentPairs.clear();
        }

        {
            OrientedBoundary<3> first("+x+y+z", "+y");
            OrientedBoundary<3> second("+x-y+z", "+y");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 9);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
            for (auto pair : StaticArray<std::pair<Size,Size>,9> {
                    { 3,  3},
                    { 4,  4},
                    {12, 12},
                    {13, 13},
                    { 5,  5},
                    {21, 21},
                    {22, 22},
                    {14, 14},
                    {23, 23}
                }) {
                CIE_TEST_CHECK(std::find(
                    coincidentPairs.begin(),
                    coincidentPairs.end(),
                    pair
                ) != coincidentPairs.end());
            }
            coincidentPairs.clear();
        }

        {
            OrientedBoundary<3> first("+x+y+z", "+z");
            OrientedBoundary<3> second("+x+y-z", "+z");

            CIE_TEST_CHECK(map.getPairCount(first, second) == 9);
            map.getPairs(first, second, std::back_inserter(coincidentPairs));
            CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
            for (auto pair : StaticArray<std::pair<Size,Size>,9> {
                    { 9,  9},
                    {10, 10},
                    {12, 12},
                    {13, 13},
                    {11, 11},
                    {14, 14},
                    {15, 15},
                    {16, 16},
                    {17, 17}
                }) {
                CIE_TEST_CHECK(std::find(
                    coincidentPairs.begin(),
                    coincidentPairs.end(),
                    pair
                ) != coincidentPairs.end());
            }
            coincidentPairs.clear();
        }
    } // 3D
}


} // namespace cie::fem
