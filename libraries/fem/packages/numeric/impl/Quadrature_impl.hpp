#ifndef CIE_FEM_NUMERIC_QUADRATURE_IMPL_HPP
#define CIE_FEM_NUMERIC_QUADRATURE_IMPL_HPP

// --- FEM Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/numeric/inc/Quadrature.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <concepts::Numeric TValue, unsigned Dimension>
template <maths::Expression TExpression>
inline void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> r_expression,
                                                   typename TExpression::Iterator it_outBegin) const
{
    DynamicArray<TValue> buffer(r_expression.size());
    this->evaluate(r_expression,
                   buffer.data(),
                   it_outBegin);
}


template <concepts::Numeric TValue, unsigned Dimension>
template <maths::Expression TExpression>
inline void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> r_expression,
                                                   typename TExpression::Iterator it_bufferBegin,
                                                   typename TExpression::Iterator it_out) const
{
    const unsigned size = r_expression.size();

    // Clear output
    std::fill(it_out,
              it_out + size,
              0);

    // Evaluate expression at quadrature points
    for (const auto& r_item : this->_nodesAndWeights) {
        // Evaluate expression into a buffer
        r_expression.evaluate(r_item.data(),
                              r_item.data() + Dimension,
                              it_bufferBegin);

        // Increment output with scaled buffer items
        const auto weight = r_item.back();
        for (unsigned i_out=0; i_out<size; ++i_out) {
            *(it_out + i_out) += weight * (*(it_bufferBegin + i_out));
        }
    } // for item in nodesAndWeights
}


} // namespace cie::fem


#endif
