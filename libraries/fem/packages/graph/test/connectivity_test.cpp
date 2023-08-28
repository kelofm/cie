// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"


namespace cie::fem {


CIE_TEST_CASE("connectivity", "[graph]")
{
    CIE_TEST_CASE_INIT("connectivity")
    using Basis = maths::Polynomial<double>;

    DynamicArray<Basis> basisFunctions {
        Basis(Basis::Coefficients {0.5, -0.5}),
        Basis(Basis::Coefficients {0.5, 0.5})
    };

    StaticArray<double,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    DynamicArray<std::pair<BoundaryID,unsigned>> connectivities;
    const auto functor = [&connectivities] (BoundaryID boundaryID, unsigned i_ansatz) {
        connectivities.emplace_back(boundaryID, i_ansatz);
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


} // namespace cie::fem
