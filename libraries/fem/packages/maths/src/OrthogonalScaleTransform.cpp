// --- FEM Includes ---
#include "packages/maths/inc/OrthogonalScaleTransform.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STL Includes ---
#include <algorithm>
#include <numeric>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
OrthogonalScaleTransformDerivative<TValue,Dimension>::OrthogonalScaleTransformDerivative() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              1);
}


template <concepts::Numeric TValue, unsigned Dimension>
OrthogonalScaleTransformDerivative<TValue,Dimension>::OrthogonalScaleTransformDerivative(Ref<const OrthogonalScaleTransform<TValue,Dimension>> rTransform) noexcept
{
    std::copy(rTransform._scales.begin(),
              rTransform._scales.end(),
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
TValue OrthogonalScaleTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator,
                                                                                 ConstIterator) const noexcept
{
    return std::accumulate(
        this->_scales.begin(),
        this->_scales.end(),
        1,
        [] (TValue left, TValue right) {return left * right;}
    );
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned OrthogonalScaleTransformDerivative<TValue,Dimension>::size() const noexcept
{
    return Dimension * Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
OrthogonalScaleTransform<TValue,Dimension>::OrthogonalScaleTransform() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              1);
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned OrthogonalScaleTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
typename OrthogonalScaleTransform<TValue,Dimension>::Derivative
OrthogonalScaleTransform<TValue,Dimension>::makeDerivative() const noexcept
{
    return Derivative(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename OrthogonalScaleTransform<TValue,Dimension>::Inverse
OrthogonalScaleTransform<TValue,Dimension>::makeInverse() const
{
    CIE_DIVISION_BY_ZERO_CHECK(!std::any(this->_scales.begin(),
                                         this->_scales.end(),
                                         [] (TValue scale) {return scale == 0;}))
    StaticArray<TValue,Dimension> inverseScales;
    std::transform(this->_scales.begin(),
                   this->_scales.end(),
                   inverseScales.begin(),
                   [](TValue scale) {return 1 / scale;});
    return OrthogonalScaleTransform(&inverseScales, (&inverseScales) + 1);
}


CIE_FEM_INSTANTIATE_TEMPLATE(OrthogonalScaleTransformDerivative)


CIE_FEM_INSTANTIATE_TEMPLATE(OrthogonalScaleTransform)


} // namespace cie::fem::maths
