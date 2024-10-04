#ifndef CIE_FEM_SCALE_TRANSLATE_TRANSFORM_IMPL_HPP
#define CIE_FEM_SCALE_TRANSLATE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy, std::transform (already included anyway)
#include <functional> // std::multiplies
#include <numeric> // std::accumulate


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ScaleTranslateTransformDerivative<TValue,Dimension>::evaluate(ConstIterator,
                                                              ConstIterator,
                                                              Iterator itOut) const noexcept
{
    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    static_assert(0 < Dimension);
    auto itScale = this->_scales.cbegin();
    *itOut++ = *itScale++;
    for (unsigned iColumn=0; iColumn<Dimension-1; ++iColumn) {
        for (unsigned iNullComponent=0; iNullComponent<Dimension; ++iNullComponent) {
            *itOut++ = static_cast<TValue>(0);
        } // for iNullComponent in range(Dimension)
        *itOut++ = *itScale++;
    } // for iColumn in range(Dimension-1)
}


template <concepts::Numeric TValue, unsigned Dimension>
TValue ScaleTranslateTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator, ConstIterator) const noexcept
{
    return std::accumulate(this->_scales.begin(),
                           this->_scales.end(),
                           static_cast<TValue>(1),
                           std::multiplies<TValue>());
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ScaleTranslateTransformDerivative<TValue,Dimension>::size() const noexcept
{
    return Dimension * Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
inline ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(TPointIt itTransformedBegin,
                                                                          [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(
        std::distance(itTransformedBegin, itTransformedEnd) == 2,
        "Expecting 2 points, but got " << std::distance(itTransformedBegin, itTransformedEnd)
    )

    const auto& rBase = *itTransformedBegin;
    const auto& rOp = *(itTransformedBegin + 1);
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        const TValue diff = rOp[iDim] - rBase[iDim];
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[iDim] = diff / static_cast<TValue>(2);
        this->_offset[iDim] = (rOp[iDim] + rBase[iDim]) / static_cast<TValue>(2);
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ScaleTranslateTransform<TValue,Dimension>::evaluate(ConstIterator itArgumentBegin,
                                                    [[maybe_unused]] ConstIterator itArgumentEnd,
                                                    Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itArgumentBegin, itArgumentEnd) == Dimension)
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        *itOut++ = itArgumentBegin[iDim] * this->_scales[iDim] + this->_offset[iDim];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              1);
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
inline TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform(TPointIt itTransformedBegin,
                                                                          [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(
        std::distance(itTransformedBegin, itTransformedEnd) == 2,
        "Expecting 2 points, but got " << std::distance(itTransformedBegin, itTransformedEnd)
    )

    const auto& rBase = *itTransformedBegin;
    const auto& rOp = *(itTransformedBegin + 1);
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        const TValue diff = (rOp[iDim] - rBase[iDim]);
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[iDim] = diff / static_cast<TValue>(2);
        this->_offset[iDim] = (rOp[iDim] + rBase[iDim]) / diff;
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
TranslateScaleTransform<TValue,Dimension>::evaluate(ConstIterator itArgumentBegin,
                                                    [[maybe_unused]] ConstIterator itArgumentEnd,
                                                    Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itArgumentBegin, itArgumentEnd) == Dimension)
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        *itOut++ = (itArgumentBegin[iDim] + this->_offset[iDim]) * this->_scales[iDim];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned TranslateScaleTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


} // namespace cie::fem::maths


#endif
