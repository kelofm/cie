#ifndef CIE_FEM_TRANSFORMED_INTEGRAND_HPP
#define CIE_FEM_TRANSFORMED_INTEGRAND_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <Expression TIntegrand, SpatialTransformDerivative TJacobian>
class TransformedIntegrand : public ExpressionTraits<typename TIntegrand::Value>
{
public:
    using typename ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename ExpressionTraits<typename TIntegrand::Value>::ConstIterator;

    using typename ExpressionTraits<typename TIntegrand::Value>::Iterator;

    using Integrand = TIntegrand;

    using Jacobian = TJacobian;

public:
    TransformedIntegrand() noexcept {}

    TransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rJacobian) noexcept
        : _integrand(std::move(rIntegrand)),
          _pJacobian(&rJacobian)
    {}

    unsigned size() const
    {return _integrand.size();}

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const
    {
        const Value scale = std::abs(_pJacobian->evaluateDeterminant(itArgumentBegin, itArgumentEnd));
        _integrand.evaluate(itArgumentBegin,
                            itArgumentEnd,
                            itOut);

        for (Value& rComponent : std::span<Value>(itOut, _integrand.size())) {
            rComponent *= scale;
        }
    }

private:
    TIntegrand _integrand;

    Ptr<const TJacobian> _pJacobian;
}; // class TransformedIntegrand


template <Expression TIntegrand, SpatialTransformDerivative TJacobian>
TransformedIntegrand<TIntegrand,TJacobian>
makeTransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         Ref<const TJacobian> rJacobian)
{
    return TransformedIntegrand(std::move(rIntegrand), rJacobian);
}


} // namespace cie::fem::maths


#endif
