#ifndef CIE_FEM_TRANSFORMED_INTEGRAND_HPP
#define CIE_FEM_TRANSFORMED_INTEGRAND_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <Expression TIntegrand, JacobianExpression TJacobian>
class TransformedIntegrand : public ExpressionTraits<typename TIntegrand::Value>
{
public:
    using typename ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename ExpressionTraits<typename TIntegrand::Value>::ConstSpan;

    using typename ExpressionTraits<typename TIntegrand::Value>::Span;

    using Integrand = TIntegrand;

    using Jacobian = TJacobian;

public:
    TransformedIntegrand() noexcept {}

    TransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rJacobian) noexcept
        : _integrand(std::move(rIntegrand)),
          _pJacobian(&rJacobian)
    {}

    unsigned size() const noexcept
    {return _integrand.size();}

    void evaluate(ConstSpan in, Span out) const
    {
        const Value scale = std::abs(_pJacobian->evaluateDeterminant(in));
        _integrand.evaluate(in, out);

        for (Value& rComponent : out) {
            rComponent *= scale;
        }
    }

private:
    TIntegrand _integrand;

    Ptr<const TJacobian> _pJacobian;
}; // class TransformedIntegrand


template <Expression TIntegrand, JacobianExpression TJacobian>
TransformedIntegrand<TIntegrand,TJacobian>
makeTransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rJacobian)
{
    return TransformedIntegrand(std::move(rIntegrand), rJacobian);
}


} // namespace cie::fem::maths


#endif
