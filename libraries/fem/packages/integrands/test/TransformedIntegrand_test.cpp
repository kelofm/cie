// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"


namespace cie::fem {


CIE_TEST_CASE("TransformedIntegrand", "[integrands]") {
    CIE_TEST_CASE_INIT("TransformedIntegrand")
    using Scalar = double;

    {
        CIE_TEST_CASE_INIT("1D")
        constexpr unsigned Dimension = 1u;
        using Function = maths::Polynomial<Scalar>;

        {
            // Integrate null, identity transform.
            const auto integrand = makeTransformedIntegrand(
                Function({0.0}),
                maths::IdentityTransform<Scalar,Dimension>().makeDerivative());

            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {
                Scalar integral;
                const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                    integrationOrder)));
                CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&integral, 1}));
                CIE_TEST_CHECK(integral == 0.0);
            }
        }

        {
            // Integrate null, scale-translate transform.
            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {

                for (unsigned segmentCount=1u; segmentCount<10u; ++segmentCount) {
                    Scalar integral = 0.0;
                    const Scalar segmentLength = 1.0 / segmentCount;

                    for (unsigned iSegment=0u; iSegment<segmentCount; ++iSegment) {
                        const std::array<std::array<Scalar,Dimension>,2> physicalCorners {
                            std::array<Scalar,Dimension> {iSegment * segmentLength},
                            std::array<Scalar,Dimension> {(iSegment + 1) * segmentLength}};

                        const auto integrand = makeTransformedIntegrand(
                            Function({0.0}),
                            maths::ScaleTranslateTransform<Scalar,Dimension>(
                                physicalCorners.begin(),
                                physicalCorners.end()).makeDerivative());

                        Scalar output;
                        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                            integrationOrder)));
                        CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&output, 1}));
                        integral += output;
                    } // for iSegment in range(segmentCount)

                    CIE_TEST_CHECK(integral == 0.0);
                } // for segmentCount in range(1, 10)
            } // for integrationOrder in range(1, 10)
        }

        {
            // Integrate constant, identity transform.
            constexpr Scalar value = 3.14;

            const auto integrand = makeTransformedIntegrand(
                Function({value}),
                maths::IdentityTransform<Scalar,Dimension>().makeDerivative());

            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {
                Scalar integral;
                const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                    integrationOrder)));
                CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&integral, 1}));
                CIE_TEST_CHECK(integral == Approx(2 * value));
            }
        }

        {
            // Integrate constant, scale-translate transform.
            constexpr Scalar value = 3.14;

            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {
                for (unsigned segmentCount=1u; segmentCount<10u; ++segmentCount) {
                    Scalar integral = 0.0;
                    const Scalar segmentLength = 1.0 / segmentCount;

                    for (unsigned iSegment=0u; iSegment<segmentCount; ++iSegment) {
                        const std::array<std::array<Scalar,Dimension>,2> physicalCorners {
                            std::array<Scalar,Dimension> {iSegment * segmentLength},
                            std::array<Scalar,Dimension> {(iSegment + 1) * segmentLength}};

                        const auto integrand = makeTransformedIntegrand(
                            Function({value}),
                            maths::ScaleTranslateTransform<Scalar,Dimension>(
                                physicalCorners.begin(),
                                physicalCorners.end()).makeInverse().makeDerivative());

                        Scalar output;
                        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                            integrationOrder)));
                        CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&output, 1}));
                        integral += output;
                    } // for iSegment in range(segmentCount)

                    CIE_TEST_CHECK(integral == Approx(value));
                } // for segmentCount in range(1, 10)
            } // for integrationOrder in range(1, 10)
        }

        {
            // Integrate linear function, identity transform.
            const Function referenceIntegral({0.0, 2.71, -3.14});
            Scalar referenceValue = 0.0;
            {
                Scalar buffer, position = 1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue += buffer;
                position = -1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue -= buffer;
            }

            const auto integrand = makeTransformedIntegrand(
                Function(referenceIntegral.makeDerivative()),
                maths::IdentityTransform<Scalar,Dimension>().makeDerivative());

            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {
                Scalar integral;
                const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                    integrationOrder)));
                CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&integral, 1}));
                CIE_TEST_CHECK(integral == Approx(referenceValue));
            }
        }

        {
            // Integrate linear function, scale-translate transform
            const Function referenceIntegral({0.0, 2.71, -3.14});
            Scalar referenceValue = 0.0;
            {
                Scalar buffer, position = 1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue += buffer;
                position = -1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue -= buffer;
            }

            for (unsigned integrationOrder=1u; integrationOrder<10u; ++integrationOrder) {
                for (unsigned segmentCount=1u; segmentCount<10u; ++segmentCount) {
                    Scalar integral = 0.0;
                    const Scalar segmentLength = 2.0 / segmentCount;

                    for (unsigned iSegment=0u; iSegment<segmentCount; ++iSegment) {
                        const std::array<std::array<Scalar,Dimension>,2> physicalCorners {
                            std::array<Scalar,Dimension> {-1.0 + iSegment * segmentLength},
                            std::array<Scalar,Dimension> {-1.0 + (iSegment + 1) * segmentLength}};

                        const auto integrand = makeTransformedIntegrand(
                            Function(referenceIntegral.makeDerivative()),
                            maths::ScaleTranslateTransform<Scalar,Dimension>(
                                physicalCorners.begin(),
                                physicalCorners.end()).makeInverse().makeDerivative());

                        Scalar output;
                        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                            integrationOrder)));
                        CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&output, 1}));
                        integral += output;
                    } // for iSegment in range(segmentCount)

                    CIE_TEST_CHECK(integral == Approx(referenceValue));
                } // for segmentCount in range(1, 10)
            } // for integrationOrder in range(1, 10)
        }

        {
            // Integrate quartic function, identity transform.
            const Function referenceIntegral({0.0, 2.71, -3.14, 1.23, -9.1});
            Scalar referenceValue = 0.0;
            {
                Scalar buffer, position = 1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue += buffer;
                position = -1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue -= buffer;
            }

            const auto integrand = makeTransformedIntegrand(
                Function(referenceIntegral.makeDerivative()),
                maths::IdentityTransform<Scalar,Dimension>().makeDerivative());

            for (unsigned integrationOrder=3u; integrationOrder<10u; ++integrationOrder) {
                Scalar integral;
                const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                    integrationOrder)));
                CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&integral, 1}));
                CIE_TEST_CHECK(integral == Approx(referenceValue));
            }
        }

        {
            // Integrate quartic function, scale-translate transform
            const Function referenceIntegral({0.0, 2.71, -3.14, 1.23, -9.1});
            Scalar referenceValue = 0.0;
            {
                Scalar buffer, position = 1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue += buffer;
                position = -1.0;
                referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
                referenceValue -= buffer;
            }

            for (unsigned integrationOrder=3u; integrationOrder<10u; ++integrationOrder) {
                for (unsigned segmentCount=1u; segmentCount<10u; ++segmentCount) {
                    Scalar integral = 0.0;
                    const Scalar segmentLength = 2.0 / segmentCount;

                    for (unsigned iSegment=0u; iSegment<segmentCount; ++iSegment) {
                        const std::array<std::array<Scalar,Dimension>,2> physicalCorners {
                            std::array<Scalar,Dimension> {-1.0 + iSegment * segmentLength},
                            std::array<Scalar,Dimension> {-1.0 + (iSegment + 1) * segmentLength}};

                        const auto integrand = makeTransformedIntegrand(
                            Function(referenceIntegral.makeDerivative()),
                            maths::ScaleTranslateTransform<Scalar,Dimension>(
                                physicalCorners.begin(),
                                physicalCorners.end()).makeInverse().makeDerivative());

                        Scalar output;
                        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(
                            integrationOrder)));
                        CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&output, 1}));
                        integral += output;
                    } // for iSegment in range(segmentCount)

                    CIE_TEST_CHECK(integral == Approx(referenceValue));
                } // for segmentCount in range(1, 10)
            } // for integrationOrder in range(1, 10)
        }
    }

    {
        CIE_TEST_CASE_INIT("1D embedded in 2D")

        struct TestExpression : public maths::ExpressionTraits<double> {
            void evaluate(
                maths::ExpressionTraits<double>::ConstSpan in,
                maths::ExpressionTraits<double>::Span out) const noexcept {
                    const double f = 1.0 + 2.0 * in.front();
                    const double g = 3.0 * std::pow(in.back(), 2) + 4.0 * std::pow(in.back(), 3);
                    out.front() = f * g;
                }

            static constexpr unsigned size() noexcept {
                return 1u;
            }
        };

        TestExpression function;

        Scalar referenceValue = 0.0;
        {
            const maths::Polynomial<double> referenceIntegral({
                0.0,
                0.0,
                0.0,
                2.0,
                std::sqrt(2.0) * 5.0,
                32.0 / 5.0});
            Scalar buffer, position = std::sqrt(2.0);
            referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
            referenceValue += buffer;
            position = -std::sqrt(2.0);
            referenceIntegral.evaluate({&position, 1}, {&buffer, 1});
            referenceValue -= buffer;
        }

        for (unsigned integrationOrder=3u; integrationOrder<10u; ++integrationOrder) {
            for (unsigned segmentCount=1u; segmentCount<10u; ++segmentCount) {
                Scalar integral = 0.0;
                const Scalar segmentLength = 2.0 / segmentCount;

                for (unsigned iSegment=0u; iSegment<segmentCount; ++iSegment) {
                    const std::array<std::array<Scalar,2>,2> physicalCorners {
                        std::array<Scalar,2> {
                            -1.0 + iSegment * segmentLength,
                            -1.0 + iSegment * segmentLength},
                        std::array<Scalar,2> {
                            -1.0 + (iSegment + 1) * segmentLength,
                            -1.0 + (iSegment + 1) * segmentLength}};

                    const auto integrand = makeTransformedIntegrand(
                        decltype(function)(function),
                        maths::AffineEmbedding<double,1u,2u>(
                            physicalCorners).makeInverse().makeDerivative());

                    Scalar output;
                    const Quadrature<Scalar,1> quadrature((GaussLegendreQuadrature<Scalar>(
                        integrationOrder)));
                    CIE_TEST_CHECK_NOTHROW(quadrature.evaluate(integrand, {&output, 1}));
                    integral += output;
                } // for iSegment in range(segmentCount)

                CIE_TEST_CHECK(integral == Approx(referenceValue));
            } // for segmentCount in range(1, 10)
        } // for integrationOrder in range(1, 10)
    }
} // CIE_TEST_CASE "TransformedIntegrand"


} // namespace cie::fem
