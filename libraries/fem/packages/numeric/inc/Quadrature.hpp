#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/QuadratureBase.hpp"
#include "packages/utilities/inc/kernel.hpp"


namespace cie::fem {


template <concepts::Numeric TValue, unsigned Dimension>
class Quadrature : public Kernel<Dimension,TValue>
{
public:
    Quadrature(Ref<const QuadratureBase<TValue>> rBase);

    Quadrature(Ref<const typename QuadratureBase<TValue>::NodeContainer> rNodes,
               Ref<const typename QuadratureBase<TValue>::WeightContainer> rWeights);

    template <maths::Expression TExpression>
    void evaluate(Ref<const TExpression> rExpression,
                  typename TExpression::Span buffer,
                  typename TExpression::Span out) const;

    template <maths::Expression TExpression>
    void evaluate(Ref<const TExpression> rExpression,
                  typename TExpression::Span out) const;

private:
    DynamicArray<StaticArray<TValue,Dimension+1>> _nodesAndWeights;
}; // class Quadrature


} // namespace cie::fem

#include "packages/numeric/impl/Quadrature_impl.hpp"
