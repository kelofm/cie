// --- Internal Includes ---
#include "packages/maths/inc/LagrangePolynomial.hpp"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"


namespace cie::fem::maths {


CIE_TEST_CASE("LagrangePolynomial", "[maths]")
{
    CIE_TEST_CASE_INIT("LagrangePolynomial")

    using Test = LagrangePolynomial<double>;
    const DynamicArray<double> nodes {0.0, 1.0/4.0, 3.0/4.0, 1.0};
    double result;

    for (Size baseIndex=0; baseIndex<nodes.size(); ++baseIndex) {
        CIE_TEST_REQUIRE_NOTHROW(Test(nodes, baseIndex));
        Test polynomial(nodes, baseIndex);

        for (Size nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex) {
            CIE_TEST_CHECK_NOTHROW(polynomial.evaluate({nodes.data() + nodeIndex, nodes.data() + nodeIndex + 1}, {&result, (&result) + 1}));
            CIE_TEST_CHECK(result == (nodeIndex == baseIndex ? Approx(1.0) : Approx(0.0).margin(1e-14)));
        }
    }
}


} // namespace cie::fem::maths
