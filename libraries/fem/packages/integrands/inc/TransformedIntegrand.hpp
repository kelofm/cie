#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
class TransformedIntegrand : public maths::ExpressionTraits<typename TIntegrand::Value>
{
public:
    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::ConstSpan;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Span;

    using Integrand = TIntegrand;

    using Jacobian = TJacobian;

public:
    TransformedIntegrand() noexcept {}

    TransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rInverseJacobian) noexcept
        : _integrand(std::move(rIntegrand)),
          _pInverseJacobian(&rInverseJacobian)
    {}

    unsigned size() const noexcept
    {return _integrand.size();}

    void evaluate(ConstSpan in, Span out) const
    {
        const Value determinant = std::abs(_pInverseJacobian->evaluateDeterminant(in));
        CIE_DIVISION_BY_ZERO_CHECK(determinant != static_cast<Value>(0));
        const Value scale = static_cast<Value>(1) / determinant;

        _integrand.evaluate(in, out);

        for (Value& rComponent : out) {
            rComponent *= scale;
        }
    }

private:
    TIntegrand _integrand;

    Ptr<const TJacobian> _pInverseJacobian;
}; // class TransformedIntegrand


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
TransformedIntegrand<TIntegrand,TJacobian>
makeTransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rJacobian)
{
    return TransformedIntegrand(std::move(rIntegrand), rJacobian);
}


} // namespace cie::fem
