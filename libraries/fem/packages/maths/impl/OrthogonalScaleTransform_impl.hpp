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
OrthogonalScaleTransformDerivative<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
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
inline OrthogonalScaleTransform<TValue,Dimension>::OrthogonalScaleTransform(TPointIt it_transformedBegin,
                                                                            TPointIt it_transformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == 1)
    Ptr<const TValue> p_begin = &(*it_transformedBegin)[0];
    std::copy(p_begin,
              p_begin + Dimension,
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
OrthogonalScaleTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                     ConstIterator it_argumentEnd,
                                                     Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)
    std::transform(it_argumentBegin,
                   it_argumentEnd,
                   this->_scales.begin(),
                   it_out,
                   [] (TValue left, TValue right) {return left * right;});
}


} // namespace cie::fem::maths


#endif
