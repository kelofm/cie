#pragma once

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/maths/inc/AffineTransform.hpp"


namespace cie::fem::maths {


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue,
          unsigned InDimension,
          unsigned OutDimension>
class AffineEmbedding {};


template <concepts::Numeric TValue,
          unsigned InDimension,
          unsigned OutDimension>
class AffineEmbeddingInverse {};


template <concepts::Numeric TValue,
          unsigned InDimension,
          unsigned OutDimension>
class AffineEmbeddingDerivative {};


template <concepts::Numeric TValue,
          unsigned InDimension,
          unsigned OutDimension>
class AffineEmbeddingInverseDerivative {};


template <concepts::Numeric TValue>
class AffineEmbedding<TValue,1u,2u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbedding)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    static constexpr unsigned InDimension = 1u;

    static constexpr unsigned OutDimension = 2u;

    using InPoint = typename Kernel<InDimension,TValue>::Point;

    using OutPoint = typename Kernel<OutDimension,TValue>::Point;

    using Derivative = AffineEmbeddingDerivative<TValue,InDimension,OutDimension>;

    using Inverse = AffineEmbeddingInverse<TValue,OutDimension,InDimension>;

    AffineEmbedding() noexcept = default;

    AffineEmbedding(std::span<const OutPoint> transformed);

    void evaluate(ConstSpan in, Span out) const;

    static constexpr unsigned size() noexcept;

    Inverse makeInverse() const;

    Derivative makeDerivative() const;

private:
    friend class AffineEmbeddingInverse<TValue,2u,1u>;

    AffineEmbedding(RightRef<AffineTransform<TValue,OutDimension>> rTransform) noexcept;

    AffineTransform<TValue,OutDimension> _transform;
}; // class AffineEmbedding


template <concepts::Numeric TValue>
class AffineEmbeddingDerivative<TValue,1u,2u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    static constexpr unsigned InDimension = 1u;

    static constexpr unsigned OutDimension = 2u;

    AffineEmbeddingDerivative() noexcept = default;

    void evaluate(ConstSpan in, Span out) const;

    static constexpr unsigned size() noexcept;

    TValue evaluateDeterminant(ConstSpan in) const;

private:
    friend class AffineEmbedding<TValue,InDimension,OutDimension>;

    AffineEmbeddingDerivative(RightRef<typename AffineTransform<TValue,OutDimension>::Derivative> rTransformDerivative) noexcept;

    typename AffineTransform<TValue,OutDimension>::Derivative _transformDerivative;
}; // class AffineEmbeddingDerivative


template <concepts::Numeric TValue>
class AffineEmbeddingInverse<TValue,2u,1u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingInverse)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    static constexpr unsigned InDimension = 2u;

    static constexpr unsigned OutDimension = 1u;

    using Inverse = AffineEmbedding<TValue,OutDimension,InDimension>;

    using Derivative = AffineEmbeddingInverseDerivative<TValue,InDimension,OutDimension>;

    AffineEmbeddingInverse() noexcept = default;

    void evaluate(ConstSpan in, Span out) const;

    static constexpr unsigned size() noexcept;

    Inverse makeInverse() const;

    Derivative makeDerivative() const;

private:
    friend class AffineEmbedding<TValue,1u,2u>;

    AffineEmbeddingInverse(RightRef<typename AffineTransform<TValue,InDimension>::Inverse> rTransform) noexcept;

    typename AffineTransform<TValue,InDimension>::Inverse _transform;
}; // class AffineEmbeddingInverse


template <concepts::Numeric TValue>
class AffineEmbeddingInverseDerivative<TValue,2u,1u> : public ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(AffineEmbeddingInverseDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    static constexpr unsigned InDimension = 2u;

    static constexpr unsigned OutDimension = 1u;

    AffineEmbeddingInverseDerivative() noexcept;

    void evaluate(ConstSpan in, Span out) const;

    static constexpr unsigned size() noexcept;

    TValue evaluateDeterminant(ConstSpan in) const;

private:
    friend class AffineEmbeddingInverse<TValue,InDimension,OutDimension>;

    AffineEmbeddingInverseDerivative(RightRef<typename AffineTransform<TValue,InDimension>::Inverse::Derivative> rTransformDerivative) noexcept;

    typename AffineTransform<TValue,InDimension>::Inverse::Derivative _transformDerivative;
}; // class AffineEmbeddingInverseDerivative


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/AffineEmbedding_impl.hpp"
