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
                                                              Iterator it_out) const noexcept
{
    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    static_assert(0 < Dimension);
    auto it_scale = this->_scales.cbegin();
    *it_out++ = *it_scale++;
    for (unsigned i_column=0; i_column<Dimension-1; ++i_column) {
        for (unsigned i_nullComponent=0; i_nullComponent<Dimension; ++i_nullComponent) {
            *it_out++ = static_cast<TValue>(0);
        } // for i_nullComponent in range(Dimension)
        *it_out++ = *it_scale++;
    } // for i_column in range(Dimension-1)
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
inline ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(TPointIt it_transformedBegin,
                                                                          [[maybe_unused]] TPointIt it_transformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == 2)

    const auto& r_base = *it_transformedBegin;
    const auto& r_op = *(it_transformedBegin + 1);
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        const TValue diff = (r_op[i_dim] - r_base[i_dim]);
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[i_dim] = diff / static_cast<TValue>(2);
        this->_offset[i_dim] = (r_op[i_dim] + r_base[i_dim]) / static_cast<TValue>(2);
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ScaleTranslateTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                    [[maybe_unused]] ConstIterator it_argumentEnd,
                                                    Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        *it_out++ = it_argumentBegin[i_dim] * this->_scales[i_dim] + this->_offset[i_dim];
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
inline TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform(TPointIt it_transformedBegin,
                                                                          [[maybe_unused]] TPointIt it_transformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == 2)

    const auto& r_base = *it_transformedBegin;
    const auto& r_op = *(it_transformedBegin + 1);
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        const TValue diff = (r_op[i_dim] - r_base[i_dim]);
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[i_dim] = diff / static_cast<TValue>(2);
        this->_offset[i_dim] = (r_op[i_dim] + r_base[i_dim]) / diff;
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
TranslateScaleTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                    [[maybe_unused]] ConstIterator it_argumentEnd,
                                                    Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        *it_out++ = (it_argumentBegin[i_dim] + this->_offset[i_dim]) * this->_scales[i_dim];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned TranslateScaleTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


} // namespace cie::fem::maths


#endif
