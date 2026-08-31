// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/utilities/inc/template_macros.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(RightRef<AffineTransform<TValue,PhysicalDimension>> rTransform) noexcept
    : _transform(std::move(rTransform))
{}


template <concepts::Numeric TValue>
typename AffineEmbedding<TValue,1u,2u>::Inverse
AffineEmbedding<TValue,1u,2u>::makeInverse() const {
    return Inverse(_transform.makeInverse());
}


template <concepts::Numeric TValue>
typename AffineEmbedding<TValue,1u,2u>::Derivative
AffineEmbedding<TValue,1u,2u>::makeDerivative() const {
    return Derivative(_transform.makeDerivative());
}


template <concepts::Numeric TValue>
AffineEmbeddingDerivative<TValue,1u,2u>::AffineEmbeddingDerivative(RightRef<typename AffineTransform<TValue,PhysicalDimension>::Derivative> rTransformDerivative) noexcept
    : _transformDerivative(std::move(rTransformDerivative))
{}


template <concepts::Numeric TValue>
AffineEmbeddingInverse<TValue,2u,1u>::AffineEmbeddingInverse(RightRef<typename AffineTransform<TValue,ParametricDimension>::Inverse> rTransform) noexcept
    : _transform(std::move(rTransform))
{}


template <concepts::Numeric TValue>
typename AffineEmbeddingInverse<TValue,2u,1u>::Derivative
AffineEmbeddingInverse<TValue,2u,1u>::makeDerivative() const {
    return Derivative(_transform.makeDerivative());
}


template <concepts::Numeric TValue>
AffineEmbeddingInverseDerivative<TValue,2u,1u>::AffineEmbeddingInverseDerivative(RightRef<typename AffineTransform<TValue,ParametricDimension>::Inverse::Derivative> rTransformDerivative) noexcept
    : _transformDerivative(std::move(rTransformDerivative))
{}


template <concepts::Numeric T>
void AffineEmbedding<T,1u,2u>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _transform);
}


template <concepts::Numeric T>
template <class TA>
void AffineEmbedding<T,1u,2u>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<AffineEmbedding> rInstance,
    TA allocator,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._transform,
            allocator);
}


template <concepts::Numeric T>
void AffineEmbeddingDerivative<T,1u,2u>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _transformDerivative);
}


template <concepts::Numeric T>
    template <class TA>
void AffineEmbeddingDerivative<T,1u,2u>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<AffineEmbeddingDerivative> rInstance,
    TA allocator,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._transformDerivative,
            allocator);
}


template <concepts::Numeric T>
void AffineEmbeddingInverse<T,2u,1u>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _transform);
}


template <concepts::Numeric T>
template <class TA>
void AffineEmbeddingInverse<T,2u,1u>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<AffineEmbeddingInverse> rInstance,
    TA allocator,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._transform,
            allocator);
}


template <concepts::Numeric T>
void AffineEmbeddingInverseDerivative<T,2u,1u>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        cie::io::BinarySerializer::serialize(
            rStream,
            _transformDerivative);
}


template <concepts::Numeric T>
template <class TA>
void AffineEmbeddingInverseDerivative<T,2u,1u>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<AffineEmbeddingInverseDerivative> rInstance,
    TA allocator,
    tags::Binary) {
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._transformDerivative,
        allocator);
}


CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbedding,1u,2u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingDerivative,1u,2u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingInverse,2u,1u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingInverseDerivative,2u,1u);


} // namespace cie::fem::maths
