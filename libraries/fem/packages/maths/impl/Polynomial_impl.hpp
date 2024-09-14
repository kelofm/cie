#ifndef CIE_FEM_MATHS_POLYNOMIAL_IMPL_HPP
#define CIE_FEM_MATHS_POLYNOMIAL_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


template <class TValue>
template <concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
Polynomial<TValue>::Polynomial(TItBegin itBegin, TItEnd itEnd)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _coefficients.reserve(std::distance(itBegin, itEnd));
    std::copy(itBegin, itEnd, std::back_inserter(_coefficients));

    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
inline void Polynomial<TValue>::evaluate(ConstIterator itArgumentBegin,
                                         ConstIterator,
                                         Iterator itResultBegin) const
{
    *itResultBegin = utils::evaluatePolynomialHorner(*itArgumentBegin,
                                                      _coefficients.begin(),
                                                      _coefficients.end());
}


} // namespace cie::fem::maths


#endif
