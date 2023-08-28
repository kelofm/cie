// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


template <class TValue>
Polynomial<TValue>::Polynomial(RightRef<Coefficients> r_coefficients) noexcept
    : _coefficients(std::move(r_coefficients))
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
            const auto it_coefficientEnd = _coefficients.end();
            TValue power = static_cast<TValue>(2);
            for (auto it_coefficient=_coefficients.begin()+2; it_coefficient!=it_coefficientEnd; ++it_coefficient, ++power)
                derivativeCoefficients.push_back(power * (*it_coefficient));
        } // if 2 < polynoialOrder
    } // if 1 < polynomialOrder

    return Polynomial(std::move(derivativeCoefficients));

    CIE_END_EXCEPTION_TRACING
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(Polynomial);


} // namespace cie::fem::maths
