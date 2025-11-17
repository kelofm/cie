#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/GraphML.hpp" // io::GraphML::Serializer, io::GraphML::Deserializer

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- STL Includes ---
#include <ostream> // std::ostream


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransform;


template <concepts::Numeric TValue, unsigned Dimension>
class TranslateScaleTransform;


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransformDerivative : public ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(ScaleTranslateTransformDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

public:
    /// @brief Identity by default.
    ScaleTranslateTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(ConstSpan in, Span out) const noexcept;

    /// @brief Evaluate the determinant of the original transform's jacobian at the provided location.
    TValue evaluateDeterminant(ConstSpan in) const noexcept;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

    template <concepts::Numeric TTNumeric, unsigned Dim>
    friend Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                                        Ref<const ScaleTranslateTransformDerivative<TTNumeric,Dim>> rObject);

private:
    friend class ScaleTranslateTransform<TValue,Dimension>;

    friend class TranslateScaleTransform<TValue,Dimension>;

    ScaleTranslateTransformDerivative(Ref<const ScaleTranslateTransform<TValue,Dimension>> rTransform) noexcept;

    ScaleTranslateTransformDerivative(Ref<const TranslateScaleTransform<TValue,Dimension>> rTransform) noexcept;

private:
    StaticArray<TValue,Dimension> _scales;
}; // class ScaleTranslateTransformDerivative



/// @brief Class representing independent scaling along orthogonal coordinate axes, followed by a translation.
/// @details Uniquely defines a mapping between axis-aligned hyperrectangles in D-dimensional space.
template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransform : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(ScaleTranslateTransform)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = ScaleTranslateTransformDerivative<TValue,Dimension>;

    using Inverse = TranslateScaleTransform<TValue,Dimension>;

public:
    /// @brief Identity transform by default.
    ScaleTranslateTransform() noexcept;

    /** @brief Construct from the transformed base @f$ [-1]^D @f$ and the vertex opposite @f$ [1]^D @f$.
     *  @details The coordinates of the input transformed vertex are identical to the
     *           scaling coefficients of the transform, after undoing the translation.
     *  @note The number of input components must match the dimension.
     */
    template <concepts::Iterator TPointIt>
    ScaleTranslateTransform(TPointIt itTransformedBegin,
                            TPointIt itTransformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components
    void evaluate(ConstSpan in, Span out) const;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the transform.
    Derivative makeDerivative() const noexcept;

    /// @brief Construct the inverse of the transform.
    Inverse makeInverse() const noexcept;

    template <concepts::Numeric TTNumeric, unsigned Dim>
    friend Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                                        Ref<const ScaleTranslateTransform<TTNumeric,Dim>> rObject);

private:
    ScaleTranslateTransform(RightRef<StaticArray<TValue,Dimension>> rScales,
                            RightRef<StaticArray<TValue,Dimension>> rOffsets) noexcept;

private:
    friend class ScaleTranslateTransformDerivative<TValue,Dimension>;

    friend class TranslateScaleTransform<TValue,Dimension>;

    StaticArray<TValue,Dimension> _scales;

    StaticArray<TValue,Dimension> _offset;
}; // class ScaleTranslateTransform



/// @brief Class representing a translation followed by an independent scaling along orthogonal coordinate axes.
/// @details Uniquely defines a mapping between axis-aligned hyperrectangles in D-dimensional space.
template <concepts::Numeric TValue, unsigned Dimension>
class TranslateScaleTransform : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(TranslateScaleTransform)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = ScaleTranslateTransformDerivative<TValue,Dimension>;

    using Inverse = ScaleTranslateTransform<TValue,Dimension>;

public:
    /// @brief Identity transform by default.
    TranslateScaleTransform() noexcept;

    /** @brief Construct from the transformed base and the vertex opposite @f$ [1]^D @f$.
     *  @details The coordinates of the input transformed vertex are identical to the
     *           scaling coefficients of the transform, after undoing the translation.
     *  @note The number of input components must match the dimension.
     */
    template <concepts::Iterator TPointIt>
    TranslateScaleTransform(TPointIt itTransformedBegin,
                            TPointIt itTransformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components
    void evaluate(ConstSpan in, Span out) const;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the transform.
    ScaleTranslateTransformDerivative<TValue,Dimension> makeDerivative() const noexcept;

    /// @brief Construct the inverse of the transform.
    ScaleTranslateTransform<TValue,Dimension> makeInverse() const noexcept;

    template <concepts::Numeric TTNumeric, unsigned Dim>
    friend Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                                        Ref<const TranslateScaleTransform<TTNumeric,Dim>> rObject);

private:
    TranslateScaleTransform(RightRef<StaticArray<TValue,Dimension>> rScales,
                            RightRef<StaticArray<TValue,Dimension>> rOffsets) noexcept;

private:
    friend class ScaleTranslateTransformDerivative<TValue,Dimension>;

    friend class ScaleTranslateTransform<TValue,Dimension>;

    StaticArray<TValue,Dimension> _scales;

    StaticArray<TValue,Dimension> _offset;
}; // class TranslateScaleTransform


/// @}


} // namespace cie::fem::maths



// --- IO --- //


namespace cie::fem::io {


template <concepts::Numeric TValue, unsigned Dimension>
struct GraphML::Serializer<maths::ScaleTranslateTransformDerivative<TValue,Dimension>>
{
    void header(Ref<XMLElement> rElement) noexcept;

    void operator()(Ref<XMLElement> rElement, Ref<const maths::ScaleTranslateTransformDerivative<TValue,Dimension>> rObject) noexcept;
}; // struct io::GraphML::Serializer<ScaleTranslateTransformDerivative<TValue,Dimension>>


template <concepts::Numeric TValue, unsigned Dimension>
struct GraphML::Serializer<maths::ScaleTranslateTransform<TValue,Dimension>>
{
    void header(Ref<XMLElement> rElement) noexcept;

    void operator()(Ref<XMLElement> rElement, Ref<const maths::ScaleTranslateTransform<TValue,Dimension>> rObject) noexcept;
}; // struct io::GraphML::Serializer<ScaleTranslateTransform<TValue,Dimension>>


template <concepts::Numeric TValue, unsigned Dimension>
struct GraphML::Deserializer<maths::ScaleTranslateTransform<TValue,Dimension>>
    : public GraphML::DeserializerBase<maths::ScaleTranslateTransform<TValue,Dimension>>
{
    using Value = maths::ScaleTranslateTransform<TValue,Dimension>;

    using GraphML::DeserializerBase<Value>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);

private:
    StaticArray<TValue,Dimension * 2> _buffer;
}; // GraphML::Deserializer<ScaleTranslateTransform>


template <concepts::Numeric TValue, unsigned Dimension>
struct GraphML::Serializer<maths::TranslateScaleTransform<TValue,Dimension>>
{
    void header(Ref<XMLElement> rElement) noexcept;

    void operator()(Ref<XMLElement> rElement, Ref<const maths::TranslateScaleTransform<TValue,Dimension>> rObject) noexcept;
}; // struct io::GraphML::Serializer<TranslateScaleTransform<TValue,Dimension>>


} // namespace cie::fem::io

#include "packages/maths/impl/ScaleTranslateTransform_impl.hpp"
