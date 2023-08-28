#ifndef CIE_MATHS_CONCEPTS_EXPRESSION_HPP
#define CIE_MATHS_CONCEPTS_EXPRESSION_HPP

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
    /// Value type to perform numerical operations on (eg: @a double).
    typename T::Value;

    /// Iterator over a contiguous array of value types.
    typename T::Iterator;

    /// Iterator over a contiguous array of const value types.
    typename T::ConstIterator;

    /// Type of the function's derivative; must also be a @a Expression.
    //typename T::Derivative;

    /// Require the evaluation through the following signature:
    /// void Expression::evaluate(ConstIterator it_begin, ConstIterator it_end, Iterator it_out)
    {
        instance.evaluate(std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::ConstIterator>(),
                          std::declval<typename T::Iterator>())
    } -> std::same_as<void>;

    /// Require a size function indicating the of scalar components returned by @a evaluate.
    {constInstance.size()} -> std::same_as<unsigned>;

    /// Require a derivative factory
    /// Note that concepts cannot be recursive so the same requirements
    /// cannot be imposed on the returned type unfortunately.
    //{instance.makeDerivative()};
}; // concept Expression


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

    virtual void evaluate(ConstIterator it_argumentBegin,
                          ConstIterator it_argumentEnd,
                          Iterator it_output) const = 0;

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

    WrappedExpression(TExpression&& r_expression) noexcept;

    WrappedExpression(const TExpression& r_expression);

    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_output) const override;

    void makeDerivative(Ref<Ptr<Derivative>> rp_derivative) const override;
}; // class WrappedExpression


} // namespace cie::fem::maths

#include "packages/maths/impl/Expression_impl.hpp"

#endif
