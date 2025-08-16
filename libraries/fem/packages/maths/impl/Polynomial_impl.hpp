#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <class TValue>
template <concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
Polynomial<TValue>::Polynomial(TItBegin itBegin, TItEnd itEnd)
{
    _coefficients.reserve(std::distance(itBegin, itEnd));
    std::copy(itBegin, itEnd, std::back_inserter(_coefficients));
}


template <class TValue>
void Polynomial<TValue>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(in.size() == 1ul)
    CIE_OUT_OF_RANGE_CHECK(out.size() == 1ul)
    out.front() = utils::evaluatePolynomialHorner(in.front(),
                                                  _coefficients.begin(),
                                                  _coefficients.end());
}


} // namespace cie::fem::maths
