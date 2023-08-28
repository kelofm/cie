// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STL Includes ---
#include <algorithm>
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
ScaleTranslateTransformDerivative<TValue,Dimension>::ScaleTranslateTransformDerivative(Ref<const ScaleTranslateTransform<TValue,Dimension>> r_transform) noexcept
{
    std::copy(r_transform._scales.begin(),
              r_transform._scales.end(),
              this->_scales.begin());
}


template <concepts::Numeric TValue, unsigned Dimension>
TValue ScaleTranslateTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator it_argumentBegin,
                                                                                ConstIterator it_argumentEnd) const noexcept
{
    return std::accumulate(
        this->_scales.begin(),
        this->_scales.end(),
        static_cast<TValue>(1),
        [] (TValue left, TValue right) {return left * right;}
    );
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ScaleTranslateTransformDerivative<TValue,Dimension>::size() const noexcept
{
    return Dimension * Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              1);
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ScaleTranslateTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransformDerivative<TValue,Dimension>
ScaleTranslateTransform<TValue,Dimension>::makeDerivative() const noexcept
{
    return ScaleTranslateTransformDerivative<TValue,Dimension>(*this);
}


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransformDerivative)


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransform)


} // namespace cie::fem::maths
