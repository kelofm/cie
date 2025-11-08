#pragma once

// --- FEM Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- STD Includes ---
#include <span>


namespace cie::fem::maths {


template <concepts::Numeric TValue>
class PolynomialView : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = PolynomialView;

    using Coefficients = DynamicArray<TValue>;

    PolynomialView() noexcept = default;

    PolynomialView(Span coefficients) noexcept;

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept;

    Derivative makeDerivative(Span buffer) const;

private:
    Span _coefficients;
}; // class PolynomialView


/// @brief @ref Expression representing a scalar polynomial.
/// @ingroup fem
template <class TValue>
class Polynomial : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = Polynomial;

    using Coefficients = DynamicArray<TValue>;

public:
    /// @brief Uninitialized by default.
    Polynomial() noexcept = default;

    Polynomial(Polynomial&&) noexcept = default;

    Polynomial(const Polynomial& rRight);

    /// @brief Construct from a container of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    Polynomial(RightRef<Coefficients> rCoefficients) noexcept;

    /// @brief Construct from a range of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    Polynomial(Span coefficients);

    Polynomial& operator=(Polynomial&&) noexcept = default;

    Polynomial& operator=(const Polynomial& rRight);

    void evaluate(ConstSpan in, Span out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the @ref Polynomial.
    /// @details The returned object is also a @ref Polynomial,
    ///          irrespective of the current polynomial order.
    Polynomial makeDerivative() const;

    std::span<const TValue> coefficients() const noexcept;

private:
    Coefficients _coefficients;

    PolynomialView<TValue> _wrapped;
}; // class Polynomial


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TValue>
struct io::GraphML::Serializer<maths::Polynomial<TValue>>
{
    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::Polynomial<TValue>> rInstance);
}; // struct GraphML::Serializer<Polynomial>


template <class TValue>
struct io::GraphML::Deserializer<maths::Polynomial<TValue>>
    : public io::GraphML::DeserializerBase<maths::Polynomial<TValue>>
{
    using io::GraphML::DeserializerBase<maths::Polynomial<TValue>>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);

private:
    typename maths::Polynomial<TValue>::Coefficients _coefficients;
}; // struct GraphML::Deserializer<Polynomial>


} // namespace cie::fem::io


#include "packages/maths/impl/Polynomial_impl.hpp"
