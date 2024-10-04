#ifndef CIE_MATHS_CONCEPTS_EXPRESSION_HPP
#define CIE_MATHS_CONCEPTS_EXPRESSION_HPP


/// @defgroup fem Finite Element Module


// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <utility> // std::declval
#include <concepts> // std::same_as


namespace cie::fem::maths {


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

    // /// @brief Type of the function's derivative; must also be a @a Expression.
    // typename T::Derivative;

    /// @details Require the evaluation through the following signature:
    ///          @code
    ///          void Expression::evaluate(ConstIterator itBegin, ConstIterator itEnd, Iterator itOut)
    ///          @endcode
    {
        instance.evaluate(std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::Iterator>())
    } -> std::same_as<void>;

    /// @brief Require a size function indicating the of scalar components returned by @a evaluate.
    {constInstance.size()} -> std::same_as<unsigned>;

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
///          1) the class must have a derivative factory computing the Jacobian of the transform.
///             The class must have an alias <tt>typename T::Derivative</tt> for the type of its derivative, which must
///             satisfy @ref cie::fem::maths::SpatialTransformDerivative "SpatialTransformDerivative". The member function
///             constructing the derivative expression must have the following signature:
///             @code{.cpp}
///             typename T::Derivative T::makeDerivative() const
///             @endcode
///
///          2) the class must have an inverse factory computing the reverse transformation
///             that must also satisfy @p SpatialTransform and must be exposed as an alias
///             <tt>typename T::Inverse</tt>. The member function constructing the inverse
///             expression musthave the following signature:
///             @code{.cpp}
///             typename T::Inverse T::makeInverse() const
///             @endcode
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


template <class TValue>
struct ExpressionTraits
{
    using Value = TValue;

    using Iterator = Ptr<TValue>;

    using ConstIterator = Ptr<const TValue>;
}; // struct Traits


template <class TValue>
struct DynamicExpression : ExpressionTraits<TValue>
{
    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

    using Derivative = DynamicExpression;

    virtual void evaluate(ConstIterator itArgumentBegin,
                          ConstIterator itArgumentEnd,
                          Iterator itOutput) const = 0;

    virtual void makeDerivative(Ref<Ptr<Derivative>> rp_derivative) const = 0;
}; // class DynamicExpression


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

private:
    WrappedExpression() = default;

    WrappedExpression(TExpression&& rExpression) noexcept;

    WrappedExpression(const TExpression& rExpression);

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOutput) const override;

    void makeDerivative(Ref<Ptr<Derivative>> rp_derivative) const override;
}; // class WrappedExpression


} // namespace cie::fem::maths

#include "packages/maths/impl/Expression_impl.hpp"

#endif
