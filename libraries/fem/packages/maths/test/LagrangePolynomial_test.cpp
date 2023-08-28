
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
        CIE_TEST_REQUIRE_NOTHROW(Test(nodes.data(), nodes.data() + nodes.size(), baseIndex));
        Test polynomial(nodes.data(), nodes.data() + nodes.size(), baseIndex);

        for (Size nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex) {
            CIE_TEST_CHECK_NOTHROW(polynomial.evaluate(makePtrTo(nodes[nodeIndex]), makePtrTo(nodes[nodeIndex])+1, makePtrTo(result)));
            CIE_TEST_CHECK(result == (nodeIndex == baseIndex ? Approx(1.0) : Approx(0.0).margin(1e-14)));
        }
    }
}


} // namespace cie::fem::maths
