#pragma once

// --- FEM Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/Serializer.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransform;


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransformDerivative : public ExpressionTraits<TValue> {
public:
    static constexpr inline unsigned ParametricDimension = Dimension;

    static constexpr inline unsigned PhysicalDimension = Dimension;

    CIE_DEFINE_CLASS_POINTERS(OrthogonalScaleTransformDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

public:
    /// @brief Identity by default.
    OrthogonalScaleTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const noexcept;

    /// @brief Evaluate the determinant of the original transform's jacobian at the provided location.
    TValue evaluateDeterminant(
        ConstSpan in,
        BufferSpan buffer) const noexcept;

    /// @brief Get the number of components written by @ref evaluate.
    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<OrthogonalScaleTransformDerivative> rInstance,
        tags::Binary tag = {});

private:
    friend class OrthogonalScaleTransform<TValue,Dimension>;

    OrthogonalScaleTransformDerivative(Ref<const OrthogonalScaleTransform<TValue,Dimension>> rTransform) noexcept;

    StaticArray<TValue,Dimension> _scales;
}; // class OrthogonalScaleTransformDerivative



/// @brief Class representing independent scaling along orthogonal coordinate axes.
/// @details Uniquely defines a mapping between axis-aligned hyperrectangles in
///          D-dimensional space, that have their base vertices at the origin.
template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransform : private ExpressionTraits<TValue> {
public:
    static constexpr inline unsigned ParametricDimension = Dimension;

    static constexpr inline unsigned PhysicalDimension = Dimension;

    CIE_DEFINE_CLASS_POINTERS(OrthogonalScaleTransform)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    using Derivative = OrthogonalScaleTransformDerivative<TValue,Dimension>;

    using Inverse = OrthogonalScaleTransform;

public:
    /// @brief Identity transform by default.
    OrthogonalScaleTransform() noexcept;

    /** @brief Construct from the transformed vertex opposite the base @f$ [1]^D @f$.
     *  @details The coordinates of the input transformed vertex are identical to the
     *           scaling coefficients of the transform.
     *  @note The number of input components must match the dimension.
     */
    template <concepts::Iterator TPointIt>
    OrthogonalScaleTransform(
        TPointIt itTransformedBegin,
        TPointIt itTransformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components
    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    /// @brief Get the number of components written by @ref evaluate.
    static constexpr unsigned size() noexcept;

    static constexpr unsigned bufferSize() noexcept;

    /// @brief Construct the derivative of the transform.
    Derivative makeDerivative() const noexcept;

    /// @brief Construct the inverse transform.
    Inverse makeInverse() const;

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<OrthogonalScaleTransform> rInstance,
        tags::Binary tag = {});

private:
    friend class OrthogonalScaleTransformDerivative<TValue,Dimension>;

    StaticArray<TValue,Dimension> _scales;
}; // class OrthogonalScaleTransform


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/OrthogonalScaleTransform_impl.hpp"
