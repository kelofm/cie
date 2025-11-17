#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::fem::maths {


template <class TLambda, concepts::Numeric TValue>
requires concepts::CallableWith<TLambda,
                                typename ExpressionTraits<TValue>::ConstSpan,
                                typename ExpressionTraits<TValue>::Span>
class LambdaExpression : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

public:
    LambdaExpression(RightRef<TLambda> rLambda,
                     unsigned size) noexcept
        : _wrapped(std::move(rLambda)),
          _size(size)
    {}

    LambdaExpression(Ref<const TLambda> rLambda,
                     unsigned size)
        : LambdaExpression(TLambda(rLambda), size)
    {}

    void evaluate(ConstSpan in, Span out) const
    {this->_wrapped(in, out);}

    unsigned size() const noexcept
    {return this->_size;}

private:
    TLambda _wrapped;

    unsigned _size;
}; // class LambdaExpression


template <class TValue, class TLambda>
LambdaExpression<TLambda,TValue> makeLambdaExpression(TLambda&& rLambda, unsigned size)
{return LambdaExpression<TLambda,TValue>(std::forward<TLambda>(rLambda), size);}


} // namespace cie::fem::maths
