#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
template <class T>
requires ct::Match<T>::template Any<TValue,PhysicalCoordinate<TValue>>
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(std::span<const std::array<T,PhysicalDimension>,2> transformed) {
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
AffineEmbedding<TValue,1u,2u>::AffineEmbedding(Ref<const std::array<std::array<T,PhysicalDimension>,2>> rTransformed)
    : AffineEmbedding(std::span<const std::array<T,PhysicalDimension>,2>(rTransformed))
{}



template <concepts::Numeric TValue>
void AffineEmbedding<TValue,1u,2u>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        CIE_OUT_OF_RANGE_CHECK(in.size() == AffineEmbedding::ParametricDimension)
        CIE_OUT_OF_RANGE_CHECK(out.size() == AffineEmbedding::PhysicalDimension)

        std::array<TValue,PhysicalDimension> augmentedIn;
        augmentedIn[0] = in[0];
        augmentedIn[1] = static_cast<TValue>(-1);
        _transform.evaluate(augmentedIn, out, buffer);
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbedding<TValue,1u,2u>::size() noexcept {
    return AffineEmbedding::PhysicalDimension;
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbedding<TValue,1u,2u>::bufferSize() noexcept {
    return AffineTransform<TValue,PhysicalDimension>::bufferSize();
}


template <concepts::Numeric TValue>
void AffineEmbeddingDerivative<TValue,1u,2u>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        std::array<TValue,PhysicalDimension> augmented;
        std::array<TValue,PhysicalDimension*PhysicalDimension> outBuffer;

        std::copy_n(
            in.data(),
            ParametricDimension,
            augmented.data());
        std::fill_n(
            augmented.data() + ParametricDimension,
            PhysicalDimension - ParametricDimension,
            static_cast<TValue>(-1));
        _transformDerivative.evaluate(in, outBuffer, buffer);
        std::copy_n(
            outBuffer.data(),
            PhysicalDimension * ParametricDimension,
            out.data());
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingDerivative<TValue,1u,2u>::size() noexcept {
    return PhysicalDimension * ParametricDimension;
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingDerivative<TValue,1u,2u>::bufferSize() noexcept {
    return AffineTransform<TValue,PhysicalDimension>::Derivative::bufferSize();
}


template <concepts::Numeric TValue>
TValue AffineEmbeddingDerivative<TValue,1u,2u>::evaluateDeterminant(
    ConstSpan in,
    BufferSpan buffer) const {
        std::array<TValue,PhysicalDimension> augmented;
        std::copy_n(
            in.data(),
            ParametricDimension,
            augmented.data());
        std::fill_n(
            augmented.data() + ParametricDimension,
            PhysicalDimension - ParametricDimension,
            static_cast<TValue>(-1));
        return _transformDerivative.evaluateDeterminant(in, buffer);
}


template <concepts::Numeric TValue>
void AffineEmbeddingInverse<TValue,2u,1u>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        assert(in.size() == AffineEmbeddingInverse::ParametricDimension);
        assert(out.size() == AffineEmbeddingInverse::PhysicalDimension);
        StaticArray<TValue,AffineEmbeddingInverse::ParametricDimension> augmented;
        _transform.evaluate(in, augmented, buffer);
        out.front() = augmented.front();
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverse<TValue,2u,1u>::size() noexcept {
    return AffineEmbeddingInverse::PhysicalDimension;
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverse<TValue,2u,1u>::bufferSize() noexcept {
    return AffineTransform<TValue,ParametricDimension>::Inverse::bufferSize();
}


template <concepts::Numeric TValue>
void AffineEmbeddingInverseDerivative<TValue,2u,1u>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        CIE_OUT_OF_RANGE_CHECK(in.size() == ParametricDimension)
        CIE_OUT_OF_RANGE_CHECK(out.size() == ParametricDimension * PhysicalDimension)

        StaticArray<TValue,ParametricDimension*ParametricDimension> outBuffer;
        _transformDerivative.evaluate(in, outBuffer, buffer);

        // Copy the top left (ParametricDimension x PhysicalDimension) submatrix of the wrapped derivative.
        for (unsigned iRow=0u; iRow<ParametricDimension; ++iRow) {
            std::copy_n(
                outBuffer.begin() + iRow * ParametricDimension,
                PhysicalDimension,
                out.begin() + iRow * PhysicalDimension);
        }
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverseDerivative<TValue,2u,1u>::size() noexcept {
    return ParametricDimension * PhysicalDimension;
}


template <concepts::Numeric TValue>
constexpr unsigned AffineEmbeddingInverseDerivative<TValue,2u,1u>::bufferSize() noexcept {
    return AffineTransform<TValue,PhysicalDimension>::Inverse::Derivative::bufferSize();
}


template <concepts::Numeric TValue>
TValue AffineEmbeddingInverseDerivative<TValue,2u,1u>::evaluateDeterminant(
    ConstSpan in,
    BufferSpan buffer) const {
        CIE_OUT_OF_RANGE_CHECK(in.size() == ParametricDimension)
        return _transformDerivative.evaluateDeterminant(in, buffer);
}


} // namespace cie::fem::maths
