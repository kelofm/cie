// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"


namespace cie::fem {


CIE_TEST_CASE("OrientedBoundary", "[graph]")
{
    CIE_TEST_CASE_INIT("OrientedBoundary")

    {
        CIE_TEST_CASE_INIT("1D")
        using Boundary = OrientedBoundary<1>;

        // ID
        CIE_TEST_CHECK(Boundary("+x", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x", "-x").id() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x", "-x").id() == BoundaryID("-x"));

        // Equality
        CIE_TEST_CHECK(Boundary("+x", "+x") == Boundary("+x", "+x"));
        CIE_TEST_CHECK(Boundary("+x", "+x") != Boundary("+x", "-x"));
        CIE_TEST_CHECK(Boundary("-x", "+x") == Boundary("-x", "+x"));
        CIE_TEST_CHECK(Boundary("-x", "+x") != Boundary("-x", "-x"));

        // Local ID
        CIE_TEST_CHECK(Boundary("+x", "+x").localID() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x", "+x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x", "-x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x", "-x").localID() == BoundaryID("+x"));

        // Size
        CIE_TEST_CHECK(Boundary("+x", "+x").size() == 1u);
        CIE_TEST_CHECK(Boundary("-x", "+x").size() == 1u);
        CIE_TEST_CHECK(Boundary("+x", "-x").size() == 1u);
        CIE_TEST_CHECK(Boundary("-x", "-x").size() == 1u);

        // Component access
        CIE_TEST_CHECK(Boundary("+x", "+x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x", "+x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x", "+x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x", "+x").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x", "-x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x", "-x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x", "-x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x", "-x").at(0) == BoundaryID("-x"));

        // Iteration
        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            for (auto axis : Boundary("+x", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }


        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            for (auto axis : Boundary("-x", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            for (auto axis : Boundary("+x", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }


        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            for (auto axis : Boundary("-x", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }
    }


    {
        CIE_TEST_CASE_INIT("2D")
        using Boundary = OrientedBoundary<2>;

        // ID
        CIE_TEST_CHECK(Boundary("+x+y", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x").id() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y").id() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y").id() == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("-x+y", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x").id() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y").id() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y").id() == BoundaryID("-y"));


        CIE_TEST_CHECK(Boundary("+x-y", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x").id() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y").id() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y").id() == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("-x-y", "+x").id() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x").id() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y").id() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y").id() == BoundaryID("-y"));

        // Equality
        CIE_TEST_CHECK(Boundary("+x+y", "+x") == Boundary("+x+y", "+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x") == Boundary("+x+y", "-x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y") == Boundary("+x+y", "+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y") == Boundary("+x+y", "-y"));

        CIE_TEST_CHECK(Boundary("-x+y", "+x") == Boundary("-x+y", "+x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x") == Boundary("-x+y", "-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y") == Boundary("-x+y", "+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y") == Boundary("-x+y", "-y"));

        CIE_TEST_CHECK(Boundary("+x-y", "+x") == Boundary("+x-y", "+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x") == Boundary("+x-y", "-x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y") == Boundary("+x-y", "+y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y") == Boundary("+x-y", "-y"));

        CIE_TEST_CHECK(Boundary("-x-y", "+x") == Boundary("-x-y", "+x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x") == Boundary("-x-y", "-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y") == Boundary("-x-y", "+y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y") == Boundary("-x-y", "-y"));

        // Local ID
        CIE_TEST_CHECK(Boundary("+x+y", "+x").localID() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y").localID() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y").localID() == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("-x+y", "+x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x").localID() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y").localID() == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y").localID() == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("+x-y", "+x").localID() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y").localID() == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y").localID() == BoundaryID("+y"));

        CIE_TEST_CHECK(Boundary("-x-y", "+x").localID() == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x").localID() == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y").localID() == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y").localID() == BoundaryID("+y"));

        // Size
        CIE_TEST_CHECK(Boundary("+x+y", "+x").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x+y", "-x").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x+y", "+y").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x+y", "-y").size() == 2u);

        CIE_TEST_CHECK(Boundary("-x+y", "+x").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x+y", "-x").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x+y", "+y").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x+y", "-y").size() == 2u);

        CIE_TEST_CHECK(Boundary("+x-y", "+x").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x-y", "-x").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x-y", "+y").size() == 2u);
        CIE_TEST_CHECK(Boundary("+x-y", "-y").size() == 2u);

        CIE_TEST_CHECK(Boundary("-x-y", "+x").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x-y", "-x").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x-y", "+y").size() == 2u);
        CIE_TEST_CHECK(Boundary("-x-y", "-y").size() == 2u);

        // Component access
        CIE_TEST_CHECK(Boundary("+x+y", "+x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+x")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y")[1] == BoundaryID("+y"));

        CIE_TEST_CHECK(Boundary("-x+y", "+x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+x")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y")[1] == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y")[1] == BoundaryID("+y"));

        CIE_TEST_CHECK(Boundary("+x-y", "+x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+x")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y")[1] == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("-x-y", "+x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+x")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y")[1] == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y")[1] == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("+x+y", "+x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+x").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-x").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "+y").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x+y", "-y").at(1) == BoundaryID("+y"));

        CIE_TEST_CHECK(Boundary("-x+y", "+x").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+x").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-x").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "+y").at(1) == BoundaryID("+y"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x+y", "-y").at(1) == BoundaryID("+y"));

        CIE_TEST_CHECK(Boundary("+x-y", "+x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+x").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-x").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "+y").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Boundary("+x-y", "-y").at(1) == BoundaryID("-y"));

        CIE_TEST_CHECK(Boundary("-x-y", "+x").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+x").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-x").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "+y").at(1) == BoundaryID("-y"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Boundary("-x-y", "-y").at(1) == BoundaryID("-y"));

        // Iteration
        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("+x+y", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("+x+y", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("+x+y", "+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("+x+y", "-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("-x+y", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("-x+y", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("-x+y", "+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("+y");
            for (auto axis : Boundary("-x+y", "-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        //

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("+x-y", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("+x-y", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("+x-y", "+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("+x-y", "-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("-x-y", "+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("-x-y", "-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("-x-y", "+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("-y");
            for (auto axis : Boundary("-x-y", "-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        // Component-wise manipulation
        {
            Boundary boundary("+x+y", "+x");
            CIE_TEST_CHECK(boundary[0] == BoundaryID("+x"));
            CIE_TEST_CHECK(boundary[1] == BoundaryID("+y"));
            CIE_TEST_CHECK(boundary.id() == BoundaryID("+x"));

            boundary[0] = BoundaryID("-y");
            CIE_TEST_CHECK(boundary == Boundary("-y+y", "+x"));
            CIE_TEST_CHECK(boundary != Boundary("+x+y", "+x"));
            CIE_TEST_CHECK(boundary[0] == BoundaryID("-y"));
            CIE_TEST_CHECK(boundary[0] != BoundaryID("+x"));
            CIE_TEST_CHECK(boundary[1] == BoundaryID("+y"));
            CIE_TEST_CHECK(boundary.id() == BoundaryID("+x"));

            boundary[1] = BoundaryID("-x");
            CIE_TEST_CHECK(boundary == Boundary("-y-x", "+x"));
            CIE_TEST_CHECK(boundary != Boundary("-y+y", "+x"));
            CIE_TEST_CHECK(boundary[0] == BoundaryID("-y"));
            CIE_TEST_CHECK(boundary[1] == BoundaryID("-x"));
            CIE_TEST_CHECK(boundary.id() == BoundaryID("+x"));
        }
    }

    {
        CIE_TEST_CASE_INIT("3D")
        using Boundary = OrientedBoundary<3>;

        // ID
        CIE_TEST_CHECK(Boundary("+x+y+z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("+x+y+z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("+x+y+z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("+x+y+z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("+x+y+z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("+x+y+z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("-x+y+z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("-x+y+z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("-x+y+z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("-x+y+z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("-x+y+z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("-x+y+z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("+x-y+z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("+x-y+z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("+x-y+z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("+x-y+z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("+x-y+z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("+x-y+z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("-x-y+z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("-x-y+z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("-x-y+z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("-x-y+z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("-x-y+z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("-x-y+z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("+x+y-z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("+x+y-z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("+x+y-z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("+x+y-z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("+x+y-z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("+x+y-z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("-x+y-z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("-x+y-z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("-x+y-z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("-x+y-z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("-x+y-z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("-x+y-z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("+x-y-z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("+x-y-z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("+x-y-z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("+x-y-z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("+x-y-z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("+x-y-z", "-z").id() == "-z");

        CIE_TEST_CHECK(Boundary("-x-y-z", "+x").id() == "+x");
        CIE_TEST_CHECK(Boundary("-x-y-z", "-x").id() == "-x");
        CIE_TEST_CHECK(Boundary("-x-y-z", "+y").id() == "+y");
        CIE_TEST_CHECK(Boundary("-x-y-z", "-y").id() == "-y");
        CIE_TEST_CHECK(Boundary("-x-y-z", "+z").id() == "+z");
        CIE_TEST_CHECK(Boundary("-x-y-z", "-z").id() == "-z");
    } // 3D
}


} // namespace cie::fem
