#ifndef CIE_FEM_SCALE_TRANSLATE_TRANSFORM_IMPL_HPP
#define CIE_FEM_SCALE_TRANSLATE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy, std::transform (already included anyway)


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ScaleTranslateTransformDerivative<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                              ConstIterator it_argumentEnd,
                                                              Iterator it_out) const noexcept
{
    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    auto it_scale = this->_scales.cbegin();
    *it_out++ = *it_scale++;
    for (unsigned i_column=1; i_column<Dimension - 1; ++i_column) {
        for (unsigned i_nullComponent=0; i_nullComponent<Dimension+1; ++i_nullComponent) {
            *it_out++ = static_cast<TValue>(0);
        } // for i_nullComponent in range(Dimension+1)
        *it_out++ = *it_scale++;
    } // for i_column in range(Dimension-2)
    *it_out = *it_scale;
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
inline ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(TPointIt it_transformedBegin,
                                                                          TPointIt it_transformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == 2)

    const auto& r_base = *it_transformedBegin;
    const auto& r_op = *(it_transformedBegin + 1);
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        this->_scales[i_dim] = (r_op[i_dim] - r_base[i_dim]) / 2;
        this->_offset[i_dim] = (r_op[i_dim] + r_base[i_dim]) / 2;
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ScaleTranslateTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                    ConstIterator it_argumentEnd,
                                                    Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        *it_out++ = *(it_argumentBegin + i_dim) * this->_scales[i_dim] + this->_offset[i_dim];
    }
}


} // namespace cie::fem::maths


#endif
