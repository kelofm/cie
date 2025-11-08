// --- FEM Includes ---
#include "packages/maths/inc/HierarchicalPolynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STD Includes ---
#include <cmath>


namespace cie::fem::maths {


template <class TValue>
HierarchicalPolynomial<TValue>::HierarchicalPolynomial(unsigned index)
    : Polynomial<TValue>()
{
    typename Polynomial<TValue>::Coefficients coefficients;

    if (index == 0u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>( 0.5);
    } else if (index == 1u) {
        coefficients.resize(2);
        coefficients.front() = static_cast<TValue>( 0.5);
        coefficients.back()  = static_cast<TValue>(-0.5);
    } else {
        coefficients.resize(index);
        std::fill_n(coefficients.data(), index, static_cast<TValue>(0));
        const TValue scale = static_cast<TValue>(1) / std::tgamma(index);
        coefficients.back() = scale;
        coefficients[index % 2 ? 1 : 0] = -scale;
    }

    static_cast<Polynomial<TValue>&>(*this) = Polynomial<TValue>(std::move(coefficients));
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(HierarchicalPolynomial);


} // namespace cie::fem::maths
