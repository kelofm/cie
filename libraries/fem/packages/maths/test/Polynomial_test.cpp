// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/maths/inc/Polynomial.hpp"


namespace cie::fem::maths {


CIE_TEST_CASE("Polynomial", "[maths]")
{
    CIE_TEST_CASE_INIT("Polynomial")

    using Test = Polynomial<double>;
    double result;

    {
        CIE_TEST_REQUIRE_NOTHROW(Test({-1, 0, 1}));
        Test polynomial({-1, 0, 1});

        const DynamicArray<std::pair<double,double>> argumentValuePairs {
            {0.0, -1.0},
            {1.0, 0.0},
            {2.0, 3.0},
            {-1.0, 0.0},
            {-2.0, 3.0}
        };
        for (const auto& [argument, reference] : argumentValuePairs) {
            CIE_TEST_CHECK_NOTHROW(polynomial.evaluate(makePtrTo(argument),
                                   makePtrTo(argument)+1,
                                   makePtrTo(result)));
            CIE_TEST_CHECK(result == Approx(reference));
        }

        CIE_TEST_REQUIRE_NOTHROW(polynomial.makeDerivative());
        const auto derivative = polynomial.makeDerivative();
        const DynamicArray<std::pair<double,double>> derivativeArgumentValuePairs {
            {0.0, 0.0},
            {1.0, 2.0},
            {2.0, 4.0},
            {-1.0, -2.0},
            {-2.0, -4.0}
        };
        for (const auto& [argument, reference] : derivativeArgumentValuePairs) {
            CIE_TEST_CHECK_NOTHROW(derivative.evaluate(makePtrTo(argument),
                                                       makePtrTo(argument) + 1,
                                                       makePtrTo(result)));
            CIE_TEST_CHECK(result == Approx(reference));
        }
    }

    // Empty coefficient list
    {
        CIE_TEST_REQUIRE_NOTHROW(Test(Test::Coefficients()));
        Test polynomial(Test::Coefficients {});
        const DynamicArray<std::pair<double,double>> argumentValuePairs {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}};
        for (const auto& [argument, reference] : argumentValuePairs) {
            CIE_TEST_CHECK_NOTHROW(polynomial.evaluate(makePtrTo(argument),
                                                       makePtrTo(argument) + 1,
                                                       makePtrTo(result)));
            CIE_TEST_CHECK(result == Approx(reference));
        }

        CIE_TEST_REQUIRE_NOTHROW(polynomial.makeDerivative());
        const auto derivative = polynomial.makeDerivative();
        for (const auto& [argument, reference] : argumentValuePairs) {
            CIE_TEST_CHECK_NOTHROW(derivative.evaluate(makePtrTo(argument),
                                                       makePtrTo(argument) + 1,
                                                       makePtrTo(result)));
            CIE_TEST_CHECK(result == Approx(reference));
        }
    }
}


} // namespace cie::fem::maths