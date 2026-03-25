#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::fem::maths {


template <class TLambda, concepts::Numeric TValue>
requires concepts::CallableWith<
    TLambda,
    typename ExpressionTraits<TValue>::ConstSpan,
    typename ExpressionTraits<TValue>::Span,
    typename ExpressionTraits<TValue>::BufferSpan>
class LambdaExpression : public ExpressionTraits<TValue> {
public:
    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::BufferSpan;

public:
    LambdaExpression(
        RightRef<TLambda> rLambda,
        unsigned size,
        unsigned bufferSize) noexcept
            : _wrapped(std::move(rLambda)),
              _size(size),
              _bufferSize(bufferSize)
    {}

    LambdaExpression(
        Ref<const TLambda> rLambda,
        unsigned size,
        unsigned bufferSize)
            : LambdaExpression(TLambda(rLambda), size, bufferSize)
    {}

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const
    {this->_wrapped(in, out, buffer);}

    unsigned size() const noexcept
    {return this->_size;}

private:
    TLambda _wrapped;

    unsigned _size, _bufferSize;
}; // class LambdaExpression


template <class TValue, class TLambda>
LambdaExpression<TLambda,TValue> makeLambdaExpression(
    TLambda&& rLambda,
    unsigned size,
    unsigned bufferSize)
{return LambdaExpression<TLambda,TValue>(std::forward<TLambda>(rLambda), size, bufferSize);}


} // namespace cie::fem::maths
