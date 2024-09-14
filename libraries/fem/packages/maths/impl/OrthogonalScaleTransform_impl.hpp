#ifndef CIE_FEM_ORTHOGONAL_SCALE_TRANSFORM_IMPL_HPP
#define CIE_FEM_ORTHOGONAL_SCALE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/OrthogonalScaleTransform.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy, std::transform (already included anyway)


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
OrthogonalScaleTransformDerivative<TValue,Dimension>::evaluate(ConstIterator,
                                                               ConstIterator,
                                                               Iterator itOut) const noexcept
{
    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    auto itScale = this->_scales.cbegin();
    *itOut++ = *itScale++;
    for (unsigned iColumn=1; iColumn<Dimension - 1; ++iColumn) {
        for (unsigned iNullComponent=0; iNullComponent<Dimension+1; ++iNullComponent) {
            *itOut++ = static_cast<TValue>(0);
        } // for iNullComponent in range(Dimension+1)
        *itOut++ = *itScale++;
    } // for iColumn in range(Dimension-2)
    *itOut = *itScale;
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
inline OrthogonalScaleTransform<TValue,Dimension>::OrthogonalScaleTransform(TPointIt itTransformedBegin,
                                                                            [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itTransformedBegin, itTransformedEnd) == 1)
    Ptr<const TValue> pBegin = &(*itTransformedBegin)[0];
    std::copy(pBegin,
              pBegin + Dimension,
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
OrthogonalScaleTransform<TValue,Dimension>::evaluate(ConstIterator itArgumentBegin,
                                                     ConstIterator itArgumentEnd,
                                                     Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itArgumentBegin, itArgumentEnd) == Dimension)
    std::transform(itArgumentBegin,
                   itArgumentEnd,
                   this->_scales.begin(),
                   itOut,
                   [] (TValue left, TValue right) {return left * right;});
}


} // namespace cie::fem::maths


#endif
