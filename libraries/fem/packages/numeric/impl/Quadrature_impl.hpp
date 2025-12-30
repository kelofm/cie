#pragma once

// --- FEM Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/numeric/inc/Quadrature.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <concepts::Numeric TValue, unsigned Dimension>
template <maths::Expression TExpression>
void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> rExpression,
                                            typename TExpression::Span out) const
{
    DynamicArray<TValue> buffer(rExpression.size());
    this->evaluate(rExpression, buffer, out);
}


template <concepts::Numeric TValue, unsigned Dimension>
template <maths::Expression TExpression>
void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> rExpression,
                                            typename TExpression::Span buffer,
                                            typename TExpression::Span out) const
{
    const unsigned size = rExpression.size();

    // Clear output
    std::fill(out.begin(), out.end(), static_cast<TValue>(0));

    // Evaluate expression at quadrature points
    for (const auto& rItem : this->_nodesAndWeights) {
        // Evaluate expression into a buffer
        rExpression.evaluate({rItem.data(), static_cast<std::size_t>(Dimension)}, buffer);

        // Increment output with scaled buffer items
        const auto weight = rItem.back();
        for (unsigned iOut=0; iOut<size; ++iOut) {
            out[iOut] += weight * buffer[iOut];
        }
    } // for item in nodesAndWeights
}


} // namespace cie::fem
