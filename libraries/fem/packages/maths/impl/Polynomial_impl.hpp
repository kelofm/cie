#ifndef CIE_FEM_MATHS_POLYNOMIAL_IMPL_HPP
#define CIE_FEM_MATHS_POLYNOMIAL_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


template <class TValue>
template <concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
Polynomial<TValue>::Polynomial(TItBegin it_begin, TItEnd it_end)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _coefficients.reserve(std::distance(it_begin, it_end));
    std::copy(it_begin, it_end, std::back_inserter(_coefficients));

    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
inline void Polynomial<TValue>::evaluate(ConstIterator it_argumentBegin,
                                         ConstIterator it_argumentEnd,
                                         Iterator it_resultBegin) const
{
    *it_resultBegin = utils::evaluatePolynomialHorner(*it_argumentBegin,
                                                      _coefficients.begin(),
                                                      _coefficients.end());
}


} // namespace cie::fem::maths


#endif
