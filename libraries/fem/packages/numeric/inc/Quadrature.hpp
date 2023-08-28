#ifndef CIE_FEM_NUMERIC_QUADRATURE_HPP
#define CIE_FEM_NUMERIC_QUADRATURE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/QuadratureBase.hpp"


namespace cie::fem {


template <concepts::Numeric TValue, unsigned Dimension>
class Quadrature
{
public:
    Quadrature(Ref<const QuadratureBase<TValue>> r_base);

    Quadrature(Ref<const typename QuadratureBase<TValue>::NodeContainer> r_nodes,
               Ref<const typename QuadratureBase<TValue>::WeightContainer> r_weights);

    template <maths::Expression TExpression>
    void evaluate(Ref<const TExpression> r_expression,
                  typename TExpression::Iterator it_bufferBegin,
                  typename TExpression::Iterator it_out) const;

    template <maths::Expression TExpression>
    void evaluate(Ref<const TExpression> r_expression,
                  typename TExpression::Iterator it_outBegin) const;

private:
    DynamicArray<StaticArray<TValue,Dimension+1>> _nodesAndWeights;
}; // class Quadrature


} // namespace cie::fem

#include "packages/numeric/impl/Quadrature_impl.hpp"

#endif
