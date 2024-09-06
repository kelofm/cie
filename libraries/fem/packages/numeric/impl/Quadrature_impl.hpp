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
inline void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> rExpression,
                                                   typename TExpression::Iterator itOutBegin) const
{
    DynamicArray<TValue> buffer(rExpression.size());
    this->evaluate(rExpression,
                   buffer.data(),
                   itOutBegin);
}


template <concepts::Numeric TValue, unsigned Dimension>
template <maths::Expression TExpression>
inline void Quadrature<TValue,Dimension>::evaluate(Ref<const TExpression> rExpression,
                                                   typename TExpression::Iterator itBufferBegin,
                                                   typename TExpression::Iterator itOut) const
{
    const unsigned size = rExpression.size();

    // Clear output
    std::fill(itOut,
              itOut + size,
              0);

    // Evaluate expression at quadrature points
    for (const auto& rItem : this->_nodesAndWeights) {
        // Evaluate expression into a buffer
        rExpression.evaluate(rItem.data(),
                              rItem.data() + Dimension,
                              itBufferBegin);

        // Increment output with scaled buffer items
        const auto weight = rItem.back();
        for (unsigned iOut=0; iOut<size; ++iOut) {
            itOut[iOut] += weight * itBufferBegin[iOut];
        }
    } // for item in nodesAndWeights
}


template <concepts::Numeric TValue, unsigned Dimension>
template <class TOutputIt>
void Quadrature<TValue,Dimension>::getIntegrationPoints(TOutputIt itOutput) const
{
    for (const auto& rItem : _nodesAndWeights) {
        typename Quadrature::Point point;
        memcpy(point.data(), rItem.data(), Dimension);
        *itOutput++ = std::move(point);
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
template <class TOutputIt>
void Quadrature<TValue,Dimension>::getIntegrationWeights(TOutputIt itOutput) const
{
    for (const auto& rItem : _nodesAndWeights) {
        *itOutput++ = rItem.back();
    }
}


} // namespace cie::fem


#endif
