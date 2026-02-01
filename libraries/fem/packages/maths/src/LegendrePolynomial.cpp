// --- FEM Includes ---
#include "packages/maths/inc/LegendrePolynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/NChooseK.hpp"

// --- STD Includes ---
#include <cmath>
#include <numeric>


namespace cie::fem::maths {


template <class TValue>
LegendrePolynomial<TValue>::LegendrePolynomial(unsigned index)
    : Polynomial<TValue>()
{
    typename Polynomial<TValue>::Coefficients coefficients(index + 1u);

    if (index == 0u) {
        coefficients.front() = 1.0;
    } else if (index == 1u) {
        coefficients.front() = static_cast<TValue>(0);
        coefficients.back()  = static_cast<TValue>(1);
    } else {
        std::fill_n(coefficients.data(), index + 1, static_cast<TValue>(0));

        StaticArray<TValue,3> factors;
        factors[1] = static_cast<TValue>(1) / std::pow(static_cast<TValue>(2),
                                                       static_cast<TValue>(index));

        for (unsigned i=0u; i<(index / 2) + 1; ++i) {
            const unsigned iMonomial = index - 2 * i;

            factors.front() = utils::NChooseK(index, i).numberOfPermutations();
            factors.back()  = index <= 2 * (index - i)
                              ? utils::NChooseK(2 * (index - i), index).numberOfPermutations()
                              : static_cast<TValue>(1);

            const TValue magnitude = std::accumulate(factors.begin(),
                                                     factors.end(),
                                                     static_cast<TValue>(1),
                                                     std::multiplies<TValue>());
            coefficients[iMonomial] = i % 2 ? -magnitude : magnitude;
        } // for iMonomial in range(index)
    }

    static_cast<Polynomial<TValue>&>(*this) = Polynomial<TValue>(std::move(coefficients));
}


template <class TValue>
ModifiedLegendrePolynomial<TValue>::ModifiedLegendrePolynomial(unsigned index)
    : Polynomial<TValue>()
{
    typename Polynomial<TValue>::Coefficients coefficients(index + 1u);

    if (index == 0u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>(-0.5);
    } else if (index == 1u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>( 0.5);
    } else {
        std::fill_n(coefficients.data(), index + 1, static_cast<TValue>(0));

        StaticArray<TValue,3> factors;
        factors[1] = static_cast<TValue>(1) / std::pow(static_cast<TValue>(2),
                                                       static_cast<TValue>(index));

        for (unsigned i=0u; i<(index / 2) + 1; ++i) {
            const unsigned iMonomial = index - 2 * i;

            factors.front() = utils::NChooseK(index, i).numberOfPermutations();
            factors.back()  = index <= 2 * (index - i)
                              ? utils::NChooseK(2 * (index - i), index).numberOfPermutations()
                              : static_cast<TValue>(1);

            const TValue magnitude = std::accumulate(factors.begin(),
                                                     factors.end(),
                                                     static_cast<TValue>(1),
                                                     std::multiplies<TValue>());
            coefficients[iMonomial] = i % 2 ? -magnitude : magnitude;
        } // for iMonomial in range(index)
    }

    static_cast<Polynomial<TValue>&>(*this) = Polynomial<TValue>(std::move(coefficients));
}


template <class TValue>
IntegratedLegendrePolynomial<TValue>::IntegratedLegendrePolynomial(unsigned index)
    : Polynomial<TValue>()
{
    typename Polynomial<TValue>::Coefficients coefficients;

    // Special cases.
    if (index == 0u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>(-0.5);
    } else if (index == 1u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>( 0.5);
    } else {
        coefficients.resize(index + 1);
        coefficients.front() = static_cast<TValue>(0); // <== temporary

        LegendrePolynomial<TValue> derivative(index - 1);
        for (unsigned iCoefficient=1u; iCoefficient<coefficients.size(); ++iCoefficient) {
            coefficients[iCoefficient] = derivative.coefficients()[iCoefficient - 1u] / iCoefficient;
        } // for iCoefficient in range(1, coefficients.size())

        PolynomialView<TValue> homogeneous(coefficients);
        const TValue in = static_cast<TValue>(1);
        TValue out;
        homogeneous.evaluate({&in, (&in) + 1}, {&out, (&out) + 1});
        coefficients.front() = -out;
    }

    static_cast<Polynomial<TValue>&>(*this) = Polynomial<TValue>(std::move(coefficients));
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(LegendrePolynomial);
CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(ModifiedLegendrePolynomial);
CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(IntegratedLegendrePolynomial);


} // namespace cie::fem::maths
