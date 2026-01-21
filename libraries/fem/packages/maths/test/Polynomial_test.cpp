// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/maths/inc/Polynomial.hpp"


namespace cie::fem::maths {


CIE_TEST_CASE("Polynomial", "[maths]")
{
    CIE_TEST_CASE_INIT("Polynomial")

    {
        CIE_TEST_CASE_INIT("dynamic")
        double result;

        {
            const Polynomial<double>::Coefficients coefficients {-1, 0, 1};
            const DynamicArray<std::pair<double,double>> argumentValuePairs {
                {0.0, -1.0},
                {1.0, 0.0},
                {2.0, 3.0},
                {-1.0, 0.0},
                {-2.0, 3.0}
            };

            #define CIE_TMP_CHECK(POLYNOMIAL)                                       \
                for (const auto& [argument, reference] : argumentValuePairs) {      \
                    CIE_TEST_CHECK_NOTHROW(POLYNOMIAL.evaluate(                     \
                        {&argument, (&argument) + 1},                               \
                        {&result, (&result) + 1}));                                 \
                    CIE_TEST_CHECK(result == Approx(reference));                    \
                }

            // Check construction.
            CIE_TEST_REQUIRE_NOTHROW(Polynomial<double>(coefficients));
            Polynomial<double> polynomial(coefficients);

            // Check evaluation.
            CIE_TMP_CHECK(polynomial)

            // Check move constructor.
            {
                auto pSwap = std::make_unique<Polynomial<double>>(coefficients);
                Polynomial<double> other = std::move(*pSwap);
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check copy constructor.
            {
                auto pSwap = std::make_unique<Polynomial<double>>(coefficients);
                Polynomial<double> other = *pSwap;
                CIE_TMP_CHECK((*pSwap))
                CIE_TMP_CHECK(other)
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check move assignment operator.
            {
                auto pSwap = std::make_unique<Polynomial<double>>(coefficients);
                Polynomial<double> other;
                other = std::move(*pSwap);
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check copy assignment operator.
            {
                auto pSwap = std::make_unique<Polynomial<double>>(coefficients);
                Polynomial<double> other;
                other = *pSwap;
                CIE_TMP_CHECK((*pSwap))
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            #undef CIE_TMP_CHECK

            // Check derivative construction.
            CIE_TEST_REQUIRE_NOTHROW(polynomial.makeDerivative());
            const auto derivative = polynomial.makeDerivative();

            // Check derivative evalutaion.
            const DynamicArray<std::pair<double,double>> derivativeArgumentValuePairs {
                { 0.0,  0.0},
                { 1.0,  2.0},
                { 2.0,  4.0},
                {-1.0, -2.0},
                {-2.0, -4.0}
            };
            for (const auto& [argument, reference] : derivativeArgumentValuePairs) {
                CIE_TEST_CHECK_NOTHROW(derivative.evaluate(
                    {&argument, (&argument) + 1},
                    {&result, (&result) + 1}));
                CIE_TEST_CHECK(result == Approx(reference));
            }
        }

        // Empty coefficient list
        {
            CIE_TEST_REQUIRE_NOTHROW(Polynomial<double>(Polynomial<double>::Coefficients()));
            Polynomial<double> polynomial(Polynomial<double>::Coefficients {});
            const DynamicArray<std::pair<double,double>> argumentValuePairs {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}};
            for (const auto& [argument, reference] : argumentValuePairs) {
                CIE_TEST_CHECK_NOTHROW(polynomial.evaluate(
                    {&argument, (&argument) + 1},
                    {&result, (&result) + 1}));
                CIE_TEST_CHECK(result == Approx(reference));
            }

            CIE_TEST_REQUIRE_NOTHROW(polynomial.makeDerivative());
            const auto derivative = polynomial.makeDerivative();
            for (const auto& [argument, reference] : argumentValuePairs) {
                CIE_TEST_CHECK_NOTHROW(derivative.evaluate(
                    {&argument, (&argument) + 1},
                    {&result, (&result) + 1}));
                CIE_TEST_CHECK(result == Approx(reference));
            }
        }
    }

    {
        CIE_TEST_CASE_INIT("static")
        double result;

        {
            const Polynomial<double,2>::Coefficients coefficients {-1, 0, 1};
            const DynamicArray<std::pair<double,double>> argumentValuePairs {
                {0.0, -1.0},
                {1.0, 0.0},
                {2.0, 3.0},
                {-1.0, 0.0},
                {-2.0, 3.0}
            };

            #define CIE_TMP_CHECK(POLYNOMIAL)                                       \
                for (const auto& [argument, reference] : argumentValuePairs) {      \
                    CIE_TEST_CHECK_NOTHROW(POLYNOMIAL.evaluate(                     \
                        {&argument, (&argument) + 1},                               \
                        {&result, (&result) + 1}));                                 \
                    CIE_TEST_CHECK(result == Approx(reference));                    \
                }

            // Check construction.
            CIE_TEST_REQUIRE_NOTHROW(Polynomial<double,2>(coefficients));
            Polynomial<double,2> polynomial(coefficients);

            // Check evaluation.
            CIE_TMP_CHECK(polynomial)

            // Check move constructor.
            {
                auto pSwap = std::make_unique<Polynomial<double,2>>(coefficients);
                Polynomial<double,2> other = std::move(*pSwap);
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check copy constructor.
            {
                auto pSwap = std::make_unique<Polynomial<double,2>>(coefficients);
                Polynomial<double,2> other = *pSwap;
                CIE_TMP_CHECK((*pSwap))
                CIE_TMP_CHECK(other)
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check move assignment operator.
            {
                auto pSwap = std::make_unique<Polynomial<double,2>>(coefficients);
                Polynomial<double,2> other;
                other = std::move(*pSwap);
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            // Check copy assignment operator.
            {
                auto pSwap = std::make_unique<Polynomial<double,2>>(coefficients);
                Polynomial<double,2> other;
                other = *pSwap;
                CIE_TMP_CHECK((*pSwap))
                pSwap.reset();
                CIE_TMP_CHECK(other)
            }

            #undef CIE_TMP_CHECK

            // Check derivative construction.
            CIE_TEST_REQUIRE_NOTHROW(polynomial.makeDerivative());
            const auto derivative = polynomial.makeDerivative();

            // Check derivative evalutaion.
            const DynamicArray<std::pair<double,double>> derivativeArgumentValuePairs {
                { 0.0,  0.0},
                { 1.0,  2.0},
                { 2.0,  4.0},
                {-1.0, -2.0},
                {-2.0, -4.0}
            };
            for (const auto& [argument, reference] : derivativeArgumentValuePairs) {
                CIE_TEST_CHECK_NOTHROW(derivative.evaluate(
                    {&argument, (&argument) + 1},
                    {&result, (&result) + 1}));
                CIE_TEST_CHECK(result == Approx(reference));
            }
        }
    }
}


} // namespace cie::fem::maths
