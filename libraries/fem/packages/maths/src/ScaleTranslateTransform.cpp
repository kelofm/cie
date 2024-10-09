// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STL Includes ---
#include <algorithm>
#include <functional>
#include <numeric>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransformDerivative<TValue,Dimension>::ScaleTranslateTransformDerivative() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              static_cast<TValue>(1));
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransformDerivative<TValue,Dimension>::ScaleTranslateTransformDerivative(Ref<const ScaleTranslateTransform<TValue,Dimension>> rTransform) noexcept
{
    std::copy(rTransform._scales.begin(),
              rTransform._scales.end(),
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransformDerivative<TValue,Dimension>::ScaleTranslateTransformDerivative(Ref<const TranslateScaleTransform<TValue,Dimension>> rTransform) noexcept
{
    std::copy(rTransform._scales.begin(),
              rTransform._scales.end(),
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(RightRef<StaticArray<TValue,Dimension>> rScales,
                                                                   RightRef<StaticArray<TValue,Dimension>> rOffset) noexcept
    : _scales(std::move(rScales)),
      _offset(std::move(rOffset))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ScaleTranslateTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
typename ScaleTranslateTransform<TValue,Dimension>::Derivative
ScaleTranslateTransform<TValue,Dimension>::makeDerivative() const noexcept
{
    return ScaleTranslateTransformDerivative<TValue,Dimension>(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename ScaleTranslateTransform<TValue,Dimension>::Inverse
ScaleTranslateTransform<TValue,Dimension>::makeInverse() const noexcept
{
    StaticArray<TValue,Dimension> scales, offset;
    std::transform(_scales.begin(),
                   _scales.end(),
                   scales.begin(),
                   [](TValue scale){return static_cast<TValue>(1) / scale;});
    std::transform(_offset.begin(),
                   _offset.end(),
                   offset.begin(),
                   [](TValue component){return -component;});
    return TranslateScaleTransform<TValue,Dimension>(std::move(scales), std::move(offset));
}


template <concepts::Numeric TValue, unsigned Dimension>
TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform(RightRef<StaticArray<TValue,Dimension>> rScales,
                                                                   RightRef<StaticArray<TValue,Dimension>> rOffset) noexcept
    : _scales(std::move(rScales)),
      _offset(std::move(rOffset))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              1);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename TranslateScaleTransform<TValue,Dimension>::Derivative
TranslateScaleTransform<TValue,Dimension>::makeDerivative() const noexcept
{
    return ScaleTranslateTransformDerivative<TValue,Dimension>(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename TranslateScaleTransform<TValue,Dimension>::Inverse
TranslateScaleTransform<TValue,Dimension>::makeInverse() const noexcept
{
    StaticArray<TValue,Dimension> scales, offset;
    std::transform(_scales.begin(),
                   _scales.end(),
                   scales.begin(),
                   [](TValue scale){return static_cast<TValue>(1) / scale;});
    std::transform(_offset.begin(),
                   _offset.end(),
                   offset.begin(),
                   [](TValue component){return -component;});
    return ScaleTranslateTransform<TValue,Dimension>(std::move(scales), std::move(offset));
}


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransformDerivative)


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransform)


CIE_FEM_INSTANTIATE_TEMPLATE(TranslateScaleTransform)


} // namespace cie::fem::maths
