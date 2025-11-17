// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/OrientedAxes.hpp"


namespace cie::fem {


CIE_TEST_CASE("OrientedAxes", "[graph]")
{
    CIE_TEST_CASE_INIT("OrientedAxes")

    {
        CIE_TEST_CASE_INIT("1D")
        using Axes = OrientedAxes<1>;

        // Equality
        CIE_TEST_CHECK(Axes("+x") == Axes("+x"));

        // Size
        CIE_TEST_CHECK(Axes("+x").size() == 1u);
        CIE_TEST_CHECK(Axes("-x").size() == 1u);

        // Component access
        CIE_TEST_CHECK(Axes("+x")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("+x").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("-x")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Axes("-x").at(0) == BoundaryID("-x"));

        // Iteration
        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            for (auto axis : Axes("+x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            for (auto axis : Axes("-x")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }
    }


    {
        CIE_TEST_CASE_INIT("2D")
        using Axes = OrientedAxes<2>;

        // Equality
        CIE_TEST_CHECK(Axes("+x+y") == Axes("+x+y"));
        CIE_TEST_CHECK(Axes("-x+y") == Axes("-x+y"));
        CIE_TEST_CHECK(Axes("+x-y") == Axes("+x-y"));
        CIE_TEST_CHECK(Axes("-x-y") == Axes("-x-y"));

        CIE_TEST_CHECK(Axes(std::string_view("+x+y")) == Axes("+x+y"));
        CIE_TEST_CHECK(Axes(std::string_view("-x+y")) == Axes("-x+y"));
        CIE_TEST_CHECK(Axes(std::string_view("+x-y")) == Axes("+x-y"));
        CIE_TEST_CHECK(Axes(std::string_view("-x-y")) == Axes("-x-y"));

        // Size
        CIE_TEST_CHECK(Axes("+x+y").size() == 2u);
        CIE_TEST_CHECK(Axes("-x+y").size() == 2u);
        CIE_TEST_CHECK(Axes("+x-y").size() == 2u);
        CIE_TEST_CHECK(Axes("-x-y").size() == 2u);

        // Component access
        CIE_TEST_CHECK(Axes("+x+y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("+x+y")[1] == BoundaryID("+y"));

        CIE_TEST_CHECK(Axes("-x+y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Axes("-x+y")[1] == BoundaryID("+y"));

        CIE_TEST_CHECK(Axes("+x-y")[0] == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("+x-y")[1] == BoundaryID("-y"));

        CIE_TEST_CHECK(Axes("-x-y")[0] == BoundaryID("-x"));
        CIE_TEST_CHECK(Axes("-x-y")[1] == BoundaryID("-y"));

        CIE_TEST_CHECK(Axes("+x+y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("+x+y").at(1) == BoundaryID("+y"));

        CIE_TEST_CHECK(Axes("-x+y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Axes("-x+y").at(1) == BoundaryID("+y"));

        CIE_TEST_CHECK(Axes("+x-y").at(0) == BoundaryID("+x"));
        CIE_TEST_CHECK(Axes("+x-y").at(1) == BoundaryID("-y"));

        CIE_TEST_CHECK(Axes("-x-y").at(0) == BoundaryID("-x"));
        CIE_TEST_CHECK(Axes("-x-y").at(1) == BoundaryID("-y"));

        // Iteration
        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("+y");
            for (auto axis : Axes("+x+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("+y");
            for (auto axis : Axes("-x+y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("+x");
            reference.emplace_back("-y");
            for (auto axis : Axes("+x-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        {
            unsigned i = 0u;
            DynamicArray<BoundaryID> reference;
            reference.emplace_back("-x");
            reference.emplace_back("-y");
            for (auto axis : Axes("-x-y")) {
                CIE_TEST_CHECK(axis == reference[i++]);
            }
        }

        // Component-wise manipulation
        {
            Axes axes("+x+y");
            CIE_TEST_CHECK(axes[0] == BoundaryID("+x"));
            CIE_TEST_CHECK(axes[1] == BoundaryID("+y"));

            axes[0] = BoundaryID("-y");
            CIE_TEST_CHECK(axes == Axes("-y+y"));
            CIE_TEST_CHECK(axes != Axes("+x+y"));
            CIE_TEST_CHECK(axes[0] == BoundaryID("-y"));
            CIE_TEST_CHECK(axes[0] != BoundaryID("+x"));
            CIE_TEST_CHECK(axes[1] == BoundaryID("+y"));

            axes[1] = BoundaryID("-x");
            CIE_TEST_CHECK(axes == Axes("-y-x"));
            CIE_TEST_CHECK(axes != Axes("-y+y"));
            CIE_TEST_CHECK(axes[0] == BoundaryID("-y"));
            CIE_TEST_CHECK(axes[1] == BoundaryID("-x"));
        }
    }
}


} // namespace cie::fem
