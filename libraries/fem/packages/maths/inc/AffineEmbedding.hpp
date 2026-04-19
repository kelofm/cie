#pragma once

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/maths/inc/AffineTransform.hpp"


namespace cie::fem::maths {


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue,
          unsigned ParametricDimension,
          unsigned PhysicalDimension>
class AffineEmbedding {};


template <concepts::Numeric TValue,
          unsigned ParametricDimension,
          unsigned PhysicalDimension>
class AffineEmbeddingInverse {};


template <concepts::Numeric TValue,
          unsigned ParametricDimension,
          unsigned PhysicalDimension>
class AffineEmbeddingDerivative {};


template <concepts::Numeric TValue,
          unsigned ParametricDimension,
          unsigned PhysicalDimension>
class AffineEmbeddingInverseDerivative {};


template <concepts::Numeric TValue>
class AffineEmbedding<TValue,1u,2u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbedding)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    static constexpr unsigned ParametricDimension = 1u;

    static constexpr unsigned PhysicalDimension = 2u;

    using InPoint = std::array<Value,ParametricDimension>;

    using OutPoint = std::array<TValue,PhysicalDimension>;

    using Derivative = AffineEmbeddingDerivative<TValue,ParametricDimension,PhysicalDimension>;

    using Inverse = AffineEmbeddingInverse<TValue,PhysicalDimension,ParametricDimension>;

    AffineEmbedding() noexcept = default;

    template <class T>
    requires ct::Match<T>::template Any<TValue,PhysicalCoordinate<TValue>>
    AffineEmbedding(std::span<const std::array<T,PhysicalDimension>,2> transformed);

    template <class T>
    requires ct::Match<T>::template Any<TValue,PhysicalCoordinate<TValue>>
    AffineEmbedding(Ref<const std::array<std::array<T,PhysicalDimension>,2>> rTransformed);

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    Inverse makeInverse() const;

    Derivative makeDerivative() const;

private:
    friend class AffineEmbeddingInverse<TValue,2u,1u>;

    AffineEmbedding(RightRef<AffineTransform<TValue,PhysicalDimension>> rTransform) noexcept;

    AffineTransform<TValue,PhysicalDimension> _transform;
}; // class AffineEmbedding


template <concepts::Numeric TValue>
class AffineEmbeddingDerivative<TValue,1u,2u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    static constexpr unsigned ParametricDimension = 1u;

    static constexpr unsigned PhysicalDimension = 2u;

    AffineEmbeddingDerivative() noexcept = default;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    TValue evaluateDeterminant(
        ConstSpan in,
        BufferSpan buffer) const;

private:
    friend class AffineEmbedding<TValue,ParametricDimension,PhysicalDimension>;

    AffineEmbeddingDerivative(RightRef<typename AffineTransform<TValue,PhysicalDimension>::Derivative> rTransformDerivative) noexcept;

    typename AffineTransform<TValue,PhysicalDimension>::Derivative _transformDerivative;
}; // class AffineEmbeddingDerivative


template <concepts::Numeric TValue>
class AffineEmbeddingInverse<TValue,2u,1u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingInverse)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    static constexpr unsigned ParametricDimension = 2u;

    static constexpr unsigned PhysicalDimension = 1u;

    using Inverse = AffineEmbedding<TValue,PhysicalDimension,ParametricDimension>;

    using Derivative = AffineEmbeddingInverseDerivative<TValue,ParametricDimension,PhysicalDimension>;

    AffineEmbeddingInverse() noexcept = default;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan) const;

    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    Inverse makeInverse() const;

    Derivative makeDerivative() const;

private:
    friend class AffineEmbedding<TValue,1u,2u>;

    AffineEmbeddingInverse(RightRef<typename AffineTransform<TValue,ParametricDimension>::Inverse> rTransform) noexcept;

    typename AffineTransform<TValue,ParametricDimension>::Inverse _transform;
}; // class AffineEmbeddingInverse


template <concepts::Numeric TValue>
class AffineEmbeddingInverseDerivative<TValue,2u,1u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingInverseDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    static constexpr unsigned ParametricDimension = 2u;

    static constexpr unsigned PhysicalDimension = 1u;

    AffineEmbeddingInverseDerivative() noexcept;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    TValue evaluateDeterminant(
        ConstSpan in,
        BufferSpan buffer) const;

private:
    friend class AffineEmbeddingInverse<TValue,ParametricDimension,PhysicalDimension>;

    AffineEmbeddingInverseDerivative(RightRef<typename AffineTransform<TValue,ParametricDimension>::Inverse::Derivative> rTransformDerivative) noexcept;

    typename AffineTransform<TValue,ParametricDimension>::Inverse::Derivative _transformDerivative;
}; // class AffineEmbeddingInverseDerivative


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/AffineEmbedding_impl.hpp"
