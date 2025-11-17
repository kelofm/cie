#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
void PolynomialView<TValue>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(in.size() == 1u)
    CIE_OUT_OF_RANGE_CHECK(out.size() == 1u)
    out.front() = utils::evaluatePolynomialHorner(in.front(),
                                                  _coefficients.begin(),
                                                  _coefficients.end());
}


template <concepts::Numeric TValue>
unsigned PolynomialView<TValue>::size() const noexcept
{
    return 1u;
}


template <class TValue>
void Polynomial<TValue>::evaluate(ConstSpan in, Span out) const
{
    _wrapped.evaluate(in, out);
}


template <class TValue>
unsigned Polynomial<TValue>::size() const noexcept
{
    return _wrapped.size();
}


} // namespace cie::fem::maths
