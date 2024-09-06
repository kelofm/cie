#ifndef CIE_FEM_LAMBDA_EXPRESSION_HPP
#define CIE_FEM_LAMBDA_EXPRESSION_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::fem::maths {


template <class TLambda, concepts::Numeric TValue>
requires concepts::CallableWith<TLambda,
                                typename ExpressionTraits<TValue>::ConstIterator,
                                typename ExpressionTraits<TValue>::ConstIterator,
                                typename ExpressionTraits<TValue>::Iterator>
class LambdaExpression : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::ConstIterator;

    using typename ExpressionTraits<TValue>::Iterator;

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

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const
    {this->_wrapped(itArgumentBegin, itArgumentEnd, itOut);}

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


#endif
