// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- STL Includes ---
#include <numeric> // std::inner_product


namespace cie::fem::maths {


template <concepts::Numeric TValue>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(std::span<const OutPoint> transformed)
{
    CIE_OUT_OF_RANGE_CHECK(transformed.size() == 2u)
    if (transformed.size() != 2u) {
        CIE_THROW(Exception, "Expecting 2 transformed points, but got " << transformed.size() << ".")
    }

    OutPoint segment;
    std::transform(transformed[0].begin(),
                   transformed[0].end(),
                   transformed[1].begin(),
                   segment.begin(),
                   [](auto begin, auto end) {
                    return end - begin;
                   });
    const auto segmentNorm = std::sqrt(std::inner_product(segment.begin(),
                                                          segment.end(),
                                                          segment.begin(),
                                                          static_cast<TValue>(0)));
    CIE_DIVISION_BY_ZERO_CHECK(segmentNorm)
    const auto scale = static_cast<TValue>(2) / segmentNorm;

    StaticArray<OutPoint,3> augmented;
    std::copy(transformed.begin(),
              transformed.end(),
              augmented.begin());

    augmented[2][0] = transformed[0][0] - scale * (transformed[1][1] - transformed[0][1]);
    augmented[2][1] = transformed[0][1] + scale * (transformed[1][0] - transformed[0][0]);

    CIE_BEGIN_EXCEPTION_TRACING
    _transform = AffineTransform<TValue,2u>(augmented);
    CIE_END_EXCEPTION_TRACING
}


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
