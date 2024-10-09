#ifndef CIE_MATHS_CONCEPTS_EXPRESSION_HPP
#define CIE_MATHS_CONCEPTS_EXPRESSION_HPP


/// @defgroup fem Finite Element Library


// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <utility> // declval
#include <concepts> // same_as
#include <span> // span
#include <memory> // shared_ptr


namespace cie::fem::maths {


/// @brief General static interface for classes doing numerics.
/// @details A class satisfying @p Expression models a multivariate vector function
///          with a narrow range of features. Interacting with the function is wholly based on
///          preallocated contiguous arrays of scalars. The convention is that the input array's size is
///          known by both @p T and its user, but the size of the output array may only be known
///          to @p T. A query function that exposes this information to the user is thus required.
///
///          Specifically, a class implementing @p Expression must have:
///
///          - An alias for the scalar type: @p Value
///
///          - An alias for an immutable iterator spanning the input array: @p ConstIterator.
///            It is expected to be a raw pointer (<tt>using ConstIterator = const Value*</tt>).
///
///          - An alias for a mutable iterator spanning the output array: @p Iterator.
///            It is expected to be a raw pointer (<tt>using Iterator = Value*</tt>).
///
///          - A member function querying the size of the output array with the following signature:
///            @code
///            unsigned T::size() const
///            @endcode
///
///          - A member function for the evaluation of the mathematical function at a specific point.
///            - @p itArgumentBegin Iterator pointing to the first component of the input array.
///            - @p itArgumentEnd Iterator pointing one past the last component of the input array.
///            - @p itOut Iterator pointing the the first component of the output array, where the
///                       function will begin writing its results.
///            @code
///            void T::evaluate(ConstIterator itArgumentBegin,
///                             ConstIterator itArgumentEnd,
///                             Iterator itOut) const
///            @endcode
/// @see cie::fem::maths::ExpressionTraits
/// @ingroup fem
template <class T>
concept Expression
= requires (T instance, const T constInstance)
{
    /// @brief Value type to perform numerical operations on (eg: @a double).
    typename T::Value;

    /// @brief Iterator over a contiguous array of value types.
    typename T::Iterator;

    /// @brief Iterator over a contiguous array of const value types.
    typename T::ConstIterator;

    // /// @brief Type of the function's derivative; must also be an @p Expression.
    // typename T::Derivative;

    /// @brief Require a size function indicating the of scalar components returned by @a evaluate.
    {constInstance.size()} -> concepts::UnsignedInteger;

    /// @details Require the evaluation through the following signature:
    ///          @code
    ///          void Expression::evaluate(ConstIterator itBegin, ConstIterator itEnd, Iterator itOut)
    ///          @endcode
    {
        instance.evaluate(std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::Iterator>())
    } -> std::same_as<void>;

    // /// @brief Require a derivative factory
    // /// @note Concepts cannot be recursive so the same requirements
    // ///       cannot be imposed on the returned type unfortunately.
    // {instance.makeDerivative()};
}; // concept Expression



/// @brief Static interface for the derivatives of spatial transformations between different spaces of identical dimensions.
/// @details On top of the requirements defined by @ref cie::fem::maths::Expression "Expression", @p SpatialTransformDerivative
///          adds 1 extra requirement. Namely, the class must be able to compute the determinant of the transformation's
///          derivative with the following signature:
///          @code{.cpp}
///          typename T::Value T::
///          @endcode
/// @see cie::fem::maths::SpatialTransform
/// @ingroup fem
template <class T>
concept SpatialTransformDerivative
= Expression<T> && requires (const T constInstance)
{
    {
        constInstance.evaluateDeterminant(typename T::ConstIterator(),
                                          typename T::ConstIterator())
    } -> std::same_as<typename T::Value>;
}; // concept SpatialTransformDerivative



/// @brief Static interface for spatial transformations between different spaces of identical dimensions.
///
/// @details On top of the requirements defined by @ref cie::fem::maths::Expression "Expression",
///          @p SpatialTransform adds 2 extra requirements:
///
///          - the class must have a derivative factory computing the Jacobian of the transform.
///            The class must have an alias <tt>typename T::Derivative</tt> for the type of its derivative, which must
///            satisfy @ref cie::fem::maths::SpatialTransformDerivative "SpatialTransformDerivative". The member function
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
= Expression<T> && requires (const T constInstance)
{
    /// @details Require a derivative factory. The derivative type need not be a @p SpatialTransform,
    ///          but it must satisfy @ref SpatialTransformDerivative that is used for computing
    ///          @ref IntegrandTransform "transformed integrals" (they require the Jacobian's determinant).
    {constInstance.makeDerivative()} -> SpatialTransformDerivative;

    /// @details Require an inverse factory. The inverse must also be a @p SpatialTransform, but this
    ///          requirement sadly cannot be encoded recursively in C++.
    {constInstance.makeInverse()} -> Expression;
}; // concept SpatialTransform


/// @brief Static interface for an @ref cie::fem::maths::Expression "Expression" that requires a buffer during its operation.
/// @details The expression does not own the memory the buffer occupies, nor does it have the
///          means of manipulating its size. Whoever uses a @p BufferedExpression must ensure
///          that the instance is provided with an adequately sized buffer before invoking
///          <tt>T::evaluate</tt>. The purpose of such ownership is to give the option of completely
///          eliminating memory allocations from expressions for performance or use on accelerators.
///
///          On top of the requirements laid out by @ref Expression, @p BufferedExpression adds
///          2 extra functions:
///
///          - A query for the minimum buffer size with the following signature:
///            @code{.cpp}
///            unsigned T::getMinBufferSize() const
///            @endcode
///
///          - A setter for the buffer with the following signature:
///            @code
///            void T::setBuffer(std::span<typename T::Value>)
///            @endcode
///
/// @note Any class satisfying @p BufferedExpression is not expected to own the memory of its
///       buffer but rely on external management, and thus may throw exceptions when
///       provided with insufficiently sized buffers to prevent segmentation faults.
/// @ingroup fem
template <class T>
concept BufferedExpression
= Expression<T> && requires(T instance, const T constInstance)
{
    /// @brief Require a query for the minimum required buffer size.
    {constInstance.getMinBufferSize()} -> concepts::UnsignedInteger;

    /// @brief Require setting a contiguous, volatile buffer.
    {instance.setBuffer(std::span<typename T::Value>())} -> std::same_as<void>;
}; // concept BufferedExpression


/// @brief Trait class exposing type aliases required by @ref Expression.
/// @ingroup fem
template <class TValue>
struct ExpressionTraits
{
    using Value = TValue;

    using Iterator = Ptr<TValue>;

    using ConstIterator = Ptr<const TValue>;
}; // struct Traits


/// @brief Base class for expressions with dynamic polymorphism.
/// @ingroup fem
template <class TValue>
struct DynamicExpression : ExpressionTraits<TValue>
{
    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

    using Derivative = DynamicExpression;

    unsigned size() const = 0;

    virtual void evaluate(ConstIterator itArgumentBegin,
                          ConstIterator itArgumentEnd,
                          Iterator itOutput) const = 0;

    virtual std::shared_ptr<Derivative> makeDerivative() const = 0;
}; // class DynamicExpression


/// @brief Wrapper class embedding the functionality of an @ref Expression with static polynmorphism into @ref DynamicExpression.
/// @ingroup fem
template <Expression TExpression>
class WrappedExpression : public DynamicExpression<typename TExpression::Value>
{
private:
    using Base = DynamicExpression<typename TExpression::Value>;

public:
    using typename Base::Iterator;

    using typename Base::ConstIterator;

    using typename Base::Derivative;

    using ExpressionType = TExpression;

public:
    WrappedExpression() = default;

    WrappedExpression(TExpression&& rExpression) noexcept;

    WrappedExpression(const TExpression& rExpression);

    unsigned size() const override;

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOutput) const override;

    std::shared_ptr<Derivative> makeDerivative() const override;

private:
    TExpression _wrapped;
}; // class WrappedExpression


} // namespace cie::fem::maths

#include "packages/maths/impl/Expression_impl.hpp"

#endif
