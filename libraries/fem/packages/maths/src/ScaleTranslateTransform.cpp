// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(
    RightRef<std::array<TValue,Dimension>> rScales,
    RightRef<std::array<TValue,Dimension>> rOffset) noexcept
        :   _scales(std::move(rScales)),
            _offset(std::move(rOffset))
{}


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
    std::array<TValue,Dimension> scales, offset;
    std::transform(
        _scales.begin(),
        _scales.end(),
        scales.begin(),
        [](TValue scale){return static_cast<TValue>(1) / scale;});
    std::transform(
        _offset.begin(),
        _offset.end(),
        offset.begin(),
        [](TValue component){return -component;});
    return TranslateScaleTransform<TValue,Dimension>(std::move(scales), std::move(offset));
}


template <concepts::Numeric TValue, unsigned Dimension>
TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform(
    RightRef<std::array<TValue,Dimension>> rScales,
    RightRef<std::array<TValue,Dimension>> rOffset) noexcept
        :   _scales(std::move(rScales)),
            _offset(std::move(rOffset))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform() noexcept {
    std::fill(
        this->_scales.begin(),
        this->_scales.end(),
        1);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename TranslateScaleTransform<TValue,Dimension>::Derivative
TranslateScaleTransform<TValue,Dimension>::makeDerivative() const noexcept {
    return ScaleTranslateTransformDerivative<TValue,Dimension>(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
typename TranslateScaleTransform<TValue,Dimension>::Inverse
TranslateScaleTransform<TValue,Dimension>::makeInverse() const noexcept {
    std::array<TValue,Dimension> scales, offset;
    std::transform(
        _scales.begin(),
        _scales.end(),
        scales.begin(),
        [](TValue scale){return static_cast<TValue>(1) / scale;});
    std::transform(
        _offset.begin(),
        _offset.end(),
        offset.begin(),
        [](TValue component){return -component;});
    return ScaleTranslateTransform<TValue,Dimension>(std::move(scales), std::move(offset));
}


template <concepts::Numeric T, unsigned D>
void ScaleTranslateTransformDerivative<T,D>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _scales.data(),
            _scales.size());
}


template <concepts::Numeric T, unsigned D>
void ScaleTranslateTransformDerivative<T,D>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<ScaleTranslateTransformDerivative> rInstance,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._scales.data(),
            rInstance._scales.size());
}


template <concepts::Numeric T, unsigned D>
void ScaleTranslateTransform<T,D>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _scales.data(),
            _scales.size());
        cie::io::BinarySerializer::serialize(
            rStream,
            _offset.data(),
            _offset.size());
}


template <concepts::Numeric T, unsigned D>
void ScaleTranslateTransform<T,D>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<ScaleTranslateTransform> rInstance,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._scales.data(),
            rInstance._scales.size());
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._offset.data(),
            rInstance._offset.size());
}


template <concepts::Numeric T, unsigned D>
void TranslateScaleTransform<T,D>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _scales.data(),
            _scales.size());
        cie::io::BinarySerializer::serialize(
            rStream,
            _offset.data(),
            _offset.size());
}


template <concepts::Numeric T, unsigned D>
void TranslateScaleTransform<T,D>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<TranslateScaleTransform> rInstance,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._scales.data(),
            rInstance._scales.size());
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._offset.data(),
            rInstance._offset.size());
}


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransformDerivative)


CIE_FEM_INSTANTIATE_TEMPLATE(ScaleTranslateTransform)


CIE_FEM_INSTANTIATE_TEMPLATE(TranslateScaleTransform)


} // namespace cie::fem::maths
