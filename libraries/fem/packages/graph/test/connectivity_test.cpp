// --- Utility Includes ---
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
        DynamicArray<std::pair<Size,Size>> coincidentPairs;

        CIE_TEST_CHECK(map.getPairCount("-x") == 1);
        map.getPairs("-x", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 1);
        CIE_TEST_CHECK(coincidentPairs.front().first == 0);
        CIE_TEST_CHECK(coincidentPairs.front().second == 1);
        coincidentPairs.clear();

        CIE_TEST_CHECK(map.getPairCount("+x") == 1);
        map.getPairs("+x", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 1);
        CIE_TEST_CHECK(coincidentPairs.front().first == 1);
        CIE_TEST_CHECK(coincidentPairs.front().second == 0);
        coincidentPairs.clear();
    } // 1D

    {
        CIE_TEST_CASE_INIT("2D")
        maths::AnsatzSpace<Basis,2> ansatzSpace(basisFunctions);
        const auto map = makeAnsatzMap(ansatzSpace,
                                       {samples.data(), samples.size()},
                                       utils::Comparison<double>(1e-10, 1e-10));
        DynamicArray<std::pair<Size,Size>> coincidentPairs;

        CIE_TEST_CHECK(map.getPairCount("-x") == 3);
        map.getPairs("-x", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(0, 1)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(3, 4)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(6, 7)
        ) != coincidentPairs.end());
        coincidentPairs.clear();

        CIE_TEST_CHECK(map.getPairCount("+x") == 3);
        map.getPairs("+x", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(1, 0)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(4, 3)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(7, 6)
        ) != coincidentPairs.end());
        coincidentPairs.clear();

        CIE_TEST_CHECK(map.getPairCount("-y") == 3);
        map.getPairs("-y", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(0, 3)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(1, 4)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(2, 5)
        ) != coincidentPairs.end());
        coincidentPairs.clear();

        CIE_TEST_CHECK(map.getPairCount("+y") == 3);
        map.getPairs("+y", std::back_inserter(coincidentPairs));
        CIE_TEST_REQUIRE(coincidentPairs.size() == 3);
        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(3, 0)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(4, 1)
        ) != coincidentPairs.end());

        CIE_TEST_CHECK(std::find(
            coincidentPairs.begin(),
            coincidentPairs.end(),
            std::pair<Size,Size>(5, 2)
        ) != coincidentPairs.end());
        coincidentPairs.clear();
    } // 2D

    {
       CIE_TEST_CASE_INIT("3D")
       maths::AnsatzSpace<Basis,3> ansatzSpace(basisFunctions);
       const auto map = makeAnsatzMap(ansatzSpace,
                                      {samples.data(), samples.size()},
                                      utils::Comparison<double>(1e-10, 1e-10));
       DynamicArray<std::pair<Size,Size>> coincidentPairs;

       CIE_TEST_CHECK(map.getPairCount("-x") == 9);
       map.getPairs("-x", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 0,  1},
               { 3,  4},
               { 9, 10},
               {12, 13},
               { 6,  7},
               {15, 16},
               {18, 19},
               {21, 22},
               {24, 25}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();

       CIE_TEST_CHECK(map.getPairCount("+x") == 9);
       map.getPairs("+x", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 1,  0},
               { 4,  3},
               {10,  9},
               {13, 12},
               { 7,  6},
               {16, 15},
               {19, 18},
               {22, 21},
               {25, 24}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();

       CIE_TEST_CHECK(map.getPairCount("-y") == 9);
       map.getPairs("-y", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 0,  3},
               { 1,  4},
               { 9, 12},
               {10, 13},
               { 2,  5},
               {18, 21},
               {19, 22},
               {11, 14},
               {20, 23}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();

       CIE_TEST_CHECK(map.getPairCount("+y") == 9);
       map.getPairs("+y", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 3,  0},
               { 4,  1},
               {12,  9},
               {13, 10},
               { 5,  2},
               {21, 18},
               {22, 19},
               {14, 11},
               {23, 20}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();

       CIE_TEST_CHECK(map.getPairCount("-z") == 9);
       map.getPairs("-z", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 0,  9},
               { 1, 10},
               { 3, 12},
               { 4, 13},
               { 2, 11},
               { 5, 14},
               { 6, 15},
               { 7, 16},
               { 8, 17}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();

       CIE_TEST_CHECK(map.getPairCount("+z") == 9);
       map.getPairs("+z", std::back_inserter(coincidentPairs));
       CIE_TEST_REQUIRE(coincidentPairs.size() == 9);
       for (auto pair : StaticArray<std::pair<Size,Size>,9> {
               { 9,  0},
               {10,  1},
               {12,  3},
               {13,  4},
               {11,  2},
               {14,  5},
               {15,  6},
               {16,  7},
               {17,  8}
           }) {
           CIE_TEST_CHECK(std::find(
               coincidentPairs.begin(),
               coincidentPairs.end(),
               pair
           ) != coincidentPairs.end());
       }
       coincidentPairs.clear();
    } // 3D
}


} // namespace cie::fem
