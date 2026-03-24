#pragma once


/// @defgroup fem Finite Element Library


// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <utility> // declval
#include <concepts> // same_as
#include <span> // span
#include <memory> // shared_ptr
#include <type_traits> // std::remove_cvref_t


namespace cie::fem::maths {


/// @brief General static interface for classes doing numerics.
/// @details A class satisfying @p Expression models a multivariate vector function
///          with a narrow range of features. Interacting with the function is wholly based on
///          preallocated contiguous arrays of scalars. The convention is that the input array's size is
///          known by both @p T and its user, but the size of the output array may only be known
///          to @p T. A query function that exposes this information to the user is thus required.
///          The @p Expression may require a preallocated chunk of contiguous memory for internal
///          calculations, whose size it must also provide.
///
///          Specifically, a class implementing @p Expression must have:
///          - An alias for the scalar type: @p Value
///          - An alias for an immutable span of the input array: @p ConstSpan.
///          - An alias for a mutable span of the output array: @p Span.
///          - An alias for a mutable span of bytes for internal use: @p BufferSpan.
///          - A member function querying the size of the output array with the following signature:
///            @code{.cpp}
///            unsigned T::size() const
///            @endcode
///            @code{.cpp}
///            unsigned T::bufferSize() const
///            @endcode
///          - A member function for the evaluation of the mathematical function at a specific point.
///            - @p input Span over the input array.
///            - @p output Span over the adequately sized output array.
///            - @p buffer Span over an adequately sized buffer.
///            @code{.cpp}
///            void T::evaluate(ConstSpan input, Span output, BufferSpan buffer) const
///            @endcode
/// @see cie::fem::maths::ExpressionTraits
/// @see cie::fem::maths::JacobianExpression
/// @see cie::fem::maths::SpatialTransform
/// @see cie::fem::maths::StaticExpression
/// @ingroup fem
template <class T>
concept Expression
= std::is_same_v<typename std::remove_cvref_t<T>::BufferSpan,std::span<std::byte>>
&& requires (std::remove_cvref_t<T> instance, const std::remove_cvref_t<T> constInstance) {
    /// @brief Value type to perform numerical operations on (eg: @p double).
    typename std::remove_cvref_t<T>::Value;

    /// @brief Span over a contiguous array of value types.
    typename std::remove_cvref_t<T>::Span;

    /// @brief Span over a contiguous array of immutable value types.
    typename std::remove_cvref_t<T>::ConstSpan;

    /// @brief Span over a contiguous array of immutable value types used as buffer for internal calculations.
    typename std::remove_cvref_t<T>::BufferSpan;

    /// @brief Require a size function indicating the number of scalar components returned by @p evaluate.
    {constInstance.size()} -> concepts::UnsignedInteger;

    /// @brief Require a function indicating the size of the buffer in bytes.
    {constInstance.bufferSize()} -> concepts::UnsignedInteger;

    /// @details Require the evaluation through the following signature:
    ///          @code
    ///          void Expression::evaluate(ConstSpan in, Span out, BufferSpan buffer)
    ///          @endcode
    {
        instance.evaluate(
            std::declval<typename std::remove_cvref_t<T>::ConstSpan>(),
            std::declval<typename std::remove_cvref_t<T>::Span>(),
            std::declval<typename std::remove_cvref_t<T>::BufferSpan>())
    } -> std::same_as<void>;
}; // concept Expression



/// @brief Static interface for the derivatives of spatial transformations between different spaces of identical dimensions.
/// @details On top of the requirements defined by @ref cie::fem::maths::Expression "Expression", @p JacobianExpression
///          adds 1 extra requirement. Namely, the class must be able to compute the determinant of the transformation's
///          derivative with the following signature:
///          @code{.cpp}
///          typename T::Value T::evaluateDeterminant(typename T::ConstSpan, typename T::BufferSpan)
///          @endcode
/// @see cie::fem::maths::Expression
/// @see cie::fem::maths::SpatialTransform
/// @ingroup fem
template <class T>
concept JacobianExpression
= Expression<T> && requires (const std::remove_cvref_t<T> constInstance) {
    {
        constInstance.evaluateDeterminant(
            typename std::remove_cvref_t<T>::ConstSpan(),
            typename std::remove_cvref_t<T>::BufferSpan())
    } -> std::same_as<typename std::remove_cvref_t<T>::Value>;
}; // concept JacobianExpression



/// @brief Static interface for spatial transformations between different spaces of identical dimensions.
///
/// @details On top of the requirements defined by @ref cie::fem::maths::Expression "Expression",
///          @p SpatialTransform adds 2 extra requirements:
///
///          - the class must have a derivative factory computing the Jacobian of the transform.
///            The class must have an alias <tt>typename T::Derivative</tt> for the type of its derivative, which must
///            satisfy @ref cie::fem::maths::JacobianExpression "JacobianExpression". The member function
///            constructing the derivative expression must have the following signature:
///            @code{.cpp}
///            typename T::Derivative T::makeDerivative() const
///            @endcode
///
///          - the class must have an inverse factory computing the reverse transformation
///            that must also satisfy @p SpatialTransform and must be exposed as an alias
///            <tt>typename T::Inverse</tt>. The member function constructing the inverse
///            expression must have the following signature:
///            @code{.cpp}
///            typename T::Inverse T::makeInverse() const
///            @endcode
///
///          Implemented spatial transformations include:
///          - @ref cie::fem::maths::IdentityTransform "IdentityTransform"
///          - @ref cie::fem::maths::OrthogonalScaleTransform "OrthogonalScaleTransform"
///          - @ref cie::fem::maths::ScaleTranslateTransform "ScaleTranslateTransform"
///          - @ref cie::fem::maths::TranslateScaleTransform "TranslateScaleTransform"
///          - @ref cie::fem::maths::AffineTransform "AffineTransform"
///          - @ref cie::fem::maths::ProjectiveTransform "ProjectiveTransform"
/// @ingroup fem
template <class T>
concept SpatialTransform
= Expression<T> && requires (const std::remove_cvref_t<T> constInstance) {
    /// @details Require a derivative factory. The derivative type need not be a @p SpatialTransform,
    ///          but it must satisfy @ref JacobianExpression that is used for computing
    ///          @ref IntegrandTransform "transformed integrals" (they require the Jacobian's determinant).
    {constInstance.makeDerivative()} -> JacobianExpression;

    /// @details Require an inverse factory. The inverse must also be a @p SpatialTransform, but this
    ///          requirement sadly cannot be encoded recursively in C++.
    {constInstance.makeInverse()} -> Expression;
}; // concept SpatialTransform


/// @brief Static interface for @ref Expression "expressions" with static memory requirements.
/// @details A static expression must know its output size and buffer size at compile time.
template <class T>
concept StaticExpression
= Expression<T> && requires () {
    {T::size()} -> std::same_as<unsigned>;
    {T::bufferSize()} -> std::same_as<unsigned>;
};


template <class T>
struct StaticExpressionSize {
    static inline constexpr unsigned size = 0u;
    static inline constexpr unsigned bufferSize = 0u;
};


template <StaticExpression T>
struct StaticExpressionSize<T> {
    static inline constexpr unsigned value = T::size();
    static inline constexpr unsigned bufferSize = T::bufferSize();
};


/// @brief Trait class exposing type aliases required by @ref Expression.
/// @ingroup fem
template <class TValue>
struct ExpressionTraits {
    using Value = TValue;

    using Span = std::span<TValue>;

    using ConstSpan = std::span<const TValue>;

    using BufferSpan = std::span<std::byte>;
}; // struct Traits


/// @brief Base class for expressions with dynamic polymorphism.
/// @ingroup fem
template <class TValue>
struct DynamicExpression : ExpressionTraits<TValue> {
    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    using Derivative = DynamicExpression;

    virtual unsigned size() const = 0;

    virtual unsigned bufferSize() const = 0;

    virtual void evaluate(
        ConstSpan input,
        Span output,
        BufferSpan buffer) = 0;

    virtual std::shared_ptr<Derivative> makeDerivative() const = 0;
}; // class DynamicExpression


/// @brief Wrapper class embedding the functionality of an @ref Expression with static polymorphism into @ref DynamicExpression.
/// @ingroup fem
template <Expression TExpression>
class WrappedExpression : public DynamicExpression<typename TExpression::Value> {
private:
    using Base = DynamicExpression<typename TExpression::Value>;

public:
    using typename Base::Span;

    using typename Base::ConstSpan;

    using typename Base::BufferSpan;

    using typename Base::Derivative;

    using ExpressionType = TExpression;

public:
    WrappedExpression() = default;

    WrappedExpression(TExpression&& rExpression) noexcept;

    WrappedExpression(const TExpression& rExpression);

    unsigned size() const override;

    void evaluate(
        ConstSpan input,
        Span output,
        BufferSpan buffer) override;

    std::shared_ptr<Derivative> makeDerivative() const override;

private:
    TExpression _wrapped;
}; // class WrappedExpression


} // namespace cie::fem::maths

#include "packages/maths/impl/Expression_impl.hpp"
