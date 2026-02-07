#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
template <class T>
requires ct::Match<T>::template Any<TValue,PhysicalCoordinate<TValue>>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(std::span<const std::array<T,OutDimension>,2> transformed) {
    CIE_OUT_OF_RANGE_CHECK(transformed.size() == 2u)
    if (transformed.size() != 2u) {
        CIE_THROW(Exception, "Expecting 2 transformed points, but got " << transformed.size() << ".")
    }

    OutPoint segment;
    std::transform(
        transformed[0].begin(),
        transformed[0].end(),
        transformed[1].begin(),
        segment.begin(),
        std::minus<TValue>());
    const auto segmentNorm = std::sqrt(std::inner_product(
        segment.begin(),
        segment.end(),
        segment.begin(),
        static_cast<TValue>(0)));
    CIE_DIVISION_BY_ZERO_CHECK(segmentNorm)
    const auto scale = static_cast<TValue>(2) / segmentNorm;

    std::array<typename AffineTransform<TValue,2u>::Point,3> augmented;
    for (std::size_t iPoint=0ul; iPoint<transformed.size(); ++iPoint) {
        std::copy_n(
            transformed[iPoint].data(),
            transformed[iPoint].size(),
            augmented[iPoint].data());
    }

    augmented.back()[0] = transformed[0][0] - scale * (transformed[1][1] - transformed[0][1]);
    augmented.back()[1] = transformed[0][1] + scale * (transformed[1][0] - transformed[0][0]);

    CIE_BEGIN_EXCEPTION_TRACING
    _transform = AffineTransform<TValue,2u>(augmented);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue>
template <class T>
requires ct::Match<T>::template Any<TValue,PhysicalCoordinate<TValue>>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(Ref<const std::array<std::array<T,OutDimension>,2>> rTransformed)
    : AffineEmbedding(std::span<const std::array<T,OutDimension>,2>(rTransformed))
{}



template <concepts::Numeric TValue>
void AffineEmbedding<TValue,1u,2u>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(in.size() == AffineEmbedding::InDimension)
    CIE_OUT_OF_RANGE_CHECK(out.size() == AffineEmbedding::OutDimension)

    std::array<TValue,OutDimension> augmentedIn;
    augmentedIn[0] = in[0];
    augmentedIn[1] = static_cast<TValue>(-1);
    _transform.evaluate(augmentedIn, out);
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbedding<TValue,1u,2u>::size() noexcept {
    return AffineEmbedding::OutDimension;
}


template <concepts::Numeric TValue>
void AffineEmbeddingDerivative<TValue,1u,2u>::evaluate(ConstSpan in, Span out) const {
    std::array<TValue,OutDimension> augmented;
    std::array<TValue,OutDimension*OutDimension> buffer;

    std::copy_n(
        in.data(),
        InDimension,
        augmented.data());
    std::fill_n(
        augmented.data() + InDimension,
        OutDimension - InDimension,
        static_cast<TValue>(-1));
    _transformDerivative.evaluate(in, buffer);
    std::copy_n(
        buffer.data(),
        OutDimension * InDimension,
        out.data());
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingDerivative<TValue,1u,2u>::size() noexcept {
    return OutDimension * InDimension;
}


template <concepts::Numeric TValue>
TValue AffineEmbeddingDerivative<TValue,1u,2u>::evaluateDeterminant(ConstSpan in) const {
    std::array<TValue,OutDimension> augmented;
    std::copy_n(
        in.data(),
        InDimension,
        augmented.data());
    std::fill_n(
        augmented.data() + InDimension,
        OutDimension - InDimension,
        static_cast<TValue>(-1));
    return _transformDerivative.evaluateDeterminant(in);
}


template <concepts::Numeric TValue>
void AffineEmbeddingInverse<TValue,2u,1u>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(in.size() == AffineEmbeddingInverse::InDimension)
    CIE_OUT_OF_RANGE_CHECK(out.size() == AffineEmbeddingInverse::OutDimension)

    StaticArray<TValue,AffineEmbeddingInverse::InDimension> augmented;
    _transform.evaluate(in, augmented);
    out[0] = augmented[0];
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverse<TValue,2u,1u>::size() noexcept {
    return AffineEmbeddingInverse::OutDimension;
}


template <concepts::Numeric TValue>
void AffineEmbeddingInverseDerivative<TValue,2u,1u>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(in.size() == InDimension)
    CIE_OUT_OF_RANGE_CHECK(out.size() == InDimension * OutDimension)

    StaticArray<TValue,InDimension*InDimension> buffer;
    _transformDerivative.evaluate(in, buffer);

    // Copy the top left (InDimension x OutDimension) submatrix of the wrapped derivative.
    for (unsigned iRow=0u; iRow<InDimension; ++iRow) {
        std::copy_n(
            buffer.begin() + iRow * InDimension,
            OutDimension,
            out.begin() + iRow * OutDimension);
    }
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverseDerivative<TValue,2u,1u>::size() noexcept {
    return InDimension * OutDimension;
}


template <concepts::Numeric TValue>
TValue AffineEmbeddingInverseDerivative<TValue,2u,1u>::evaluateDeterminant(ConstSpan in) const {
    CIE_OUT_OF_RANGE_CHECK(in.size() == InDimension)
    return _transformDerivative.evaluateDeterminant(in);
}


} // namespace cie::fem::maths
