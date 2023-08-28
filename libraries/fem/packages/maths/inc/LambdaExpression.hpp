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
    LambdaExpression(RightRef<TLambda> r_lambda,
                     unsigned size) noexcept
        : _wrapped(std::move(r_lambda)),
          _size(size)
    {}

    LambdaExpression(Ref<const TLambda> r_lambda,
                     unsigned size)
        : LambdaExpression(TLambda(r_lambda), size)
    {}

    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const
    {this->_wrapped(it_argumentBegin, it_argumentEnd, it_out);}

    unsigned size() const noexcept
    {return this->_size;}

private:
    TLambda _wrapped;

    unsigned _size;
}; // class LambdaExpression


template <class TValue, class TLambda>
LambdaExpression<TLambda,TValue> makeLambdaExpression(TLambda&& r_lambda, unsigned size)
{return LambdaExpression<TLambda,TValue>(std::forward<TLambda>(r_lambda), size);}


} // namespace cie::fem::maths


#endif
