// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STL Includes ---
#include <numeric> // std::inner_product


namespace cie::fem::maths {


template <concepts::Numeric TValue>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(RightRef<AffineTransform<TValue,OutDimension>> rTransform) noexcept
    : _transform(std::move(rTransform))
{
}


template <concepts::Numeric TValue>
typename AffineEmbedding<TValue,1u,2u>::Inverse
AffineEmbedding<TValue,1u,2u>::makeInverse() const
{
    return Inverse(_transform.makeInverse());
}


template <concepts::Numeric TValue>
typename AffineEmbedding<TValue,1u,2u>::Derivative
AffineEmbedding<TValue,1u,2u>::makeDerivative() const
{
    return Derivative(_transform.makeDerivative());
}


template <concepts::Numeric TValue>
AffineEmbeddingDerivative<TValue,1u,2u>::AffineEmbeddingDerivative(RightRef<typename AffineTransform<TValue,OutDimension>::Derivative> rTransformDerivative) noexcept
    : _transformDerivative(std::move(rTransformDerivative))
{
}


template <concepts::Numeric TValue>
AffineEmbeddingInverse<TValue,2u,1u>::AffineEmbeddingInverse(RightRef<typename AffineTransform<TValue,InDimension>::Inverse> rTransform) noexcept
    : _transform(std::move(rTransform))
{
}


template <concepts::Numeric TValue>
typename AffineEmbeddingInverse<TValue,2u,1u>::Derivative
AffineEmbeddingInverse<TValue,2u,1u>::makeDerivative() const
{
    return Derivative(_transform.makeDerivative());
}


template <concepts::Numeric TValue>
AffineEmbeddingInverseDerivative<TValue,2u,1u>::AffineEmbeddingInverseDerivative(RightRef<typename AffineTransform<TValue,InDimension>::Inverse::Derivative> rTransformDerivative) noexcept
    : _transformDerivative(std::move(rTransformDerivative))
{
}


CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbedding,1u,2u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingDerivative,1u,2u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingInverse,2u,1u);

CIE_FEM_INSTANTIATE_MIXED_TEMPLATE(AffineEmbeddingInverseDerivative,2u,1u);


} // namespace cie::fem::maths
