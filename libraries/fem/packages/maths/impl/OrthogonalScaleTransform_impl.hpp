#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/OrthogonalScaleTransform.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy, std::transform (already included anyway)


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
void OrthogonalScaleTransformDerivative<TValue,Dimension>::evaluate(ConstSpan, [[maybe_unused]] Span out) const noexcept
{
    if constexpr (!Dimension) return;

    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    auto itOut = out.data();
    for (unsigned i=0u; i<Dimension; ++i) {
        for (unsigned j=0u; j<Dimension; ++j) {
            *itOut++ = i == j ? this->_scales[i] : static_cast<TValue>(0);
        }
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
OrthogonalScaleTransform<TValue,Dimension>::OrthogonalScaleTransform(TPointIt itTransformedBegin,
                                                                     [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itTransformedBegin, itTransformedEnd) == 1)
    Ptr<const TValue> pBegin = &(*itTransformedBegin)[0];
    std::copy(pBegin,
              pBegin + Dimension,
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
void OrthogonalScaleTransform<TValue,Dimension>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(in.size() == Dimension)
    std::transform(in.begin(),
                   in.end(),
                   this->_scales.begin(),
                   out.begin(),
                   [] (TValue left, TValue right) {return left * right;});
}


} // namespace cie::fem::maths
