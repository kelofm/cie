#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <Expression TExpression>
WrappedExpression<TExpression>::WrappedExpression(TExpression&& rExpression) noexcept
    : _wrapped(std::move(rExpression))
{}


template <Expression TExpression>
WrappedExpression<TExpression>::WrappedExpression(const TExpression& rExpression)
    : _wrapped(rExpression)
{
}


template <Expression TExpression>
unsigned WrappedExpression<TExpression>::size() const {
    return _wrapped.size();
}


template <Expression TExpression>
unsigned WrappedExpression<TExpression>::bufferSize() const {
    return _wrapped.bufferSize();
}


template <Expression TExpression>
void WrappedExpression<TExpression>::evaluate(
    ConstSpan input,
    Span output,
    BufferSpan buffer) {
        _wrapped.evaluate(
            input,
            output,
            buffer);
}


template <Expression TExpression>
std::shared_ptr<typename WrappedExpression<TExpression>::Derivative>
WrappedExpression<TExpression>::makeDerivative() const {
    return std::shared_ptr<Derivative>(new WrappedExpression<typename TExpression::Derivative>(_wrapped.makeDerivative()));
}


} // namespace cie::fem::maths
