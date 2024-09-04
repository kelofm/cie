// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


template <class TValue>
Polynomial<TValue>::Polynomial(RightRef<Coefficients> rCoefficients) noexcept
    : _coefficients(std::move(rCoefficients))
{
}


template <class TValue>
Polynomial<TValue> Polynomial<TValue>::makeDerivative() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto polynomialOrder = _coefficients.size();
    Coefficients derivativeCoefficients;

    if (1 < polynomialOrder) [[likely]] {
        // Push first coefficient (no multiplication required)
        derivativeCoefficients.push_back(_coefficients[1]);
        if (2 < polynomialOrder) {
            derivativeCoefficients.reserve(polynomialOrder - 1u);
            const auto itCoefficientEnd = _coefficients.end();
            TValue power = static_cast<TValue>(2);
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power)
                derivativeCoefficients.push_back(power * (*itCoefficient));
        } // if 2 < polynoialOrder
    } // if 1 < polynomialOrder

    return Polynomial(std::move(derivativeCoefficients));

    CIE_END_EXCEPTION_TRACING
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(Polynomial);


} // namespace cie::fem::maths
