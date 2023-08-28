// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"


namespace cie::fem {


CIE_TEST_CASE("BoundaryID", "[graph]")
{
    CIE_TEST_CASE_INIT("BoundaryID")

    {
        BoundaryID id;

        for (unsigned i_dim=0; i_dim<4; ++i_dim) {
            CIE_TEST_CHECK(id.getDimension() == i_dim);
            CIE_TEST_CHECK(id.getDirection() == 0);
            CIE_TEST_CHECK_NOTHROW(++id);
            CIE_TEST_CHECK(id.getDimension() == i_dim);
            CIE_TEST_CHECK(id.getDirection() == 1);
            CIE_TEST_CHECK_NOTHROW(++id);
        }

        CIE_TEST_CHECK(id.getDimension() == 4);
        CIE_TEST_CHECK(id.getDirection() == 0);
    }

    {
        CIE_TEST_CHECK(BoundaryID(5, 0).getDimension() == 5);
        CIE_TEST_CHECK(BoundaryID(5, 0).getDirection() == 0);
        CIE_TEST_CHECK(BoundaryID(5, 1).getDimension() == 5);
        CIE_TEST_CHECK(BoundaryID(5, 1).getDirection() == 1);
    }
}


} // namespace cie::fem
