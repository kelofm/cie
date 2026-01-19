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


/// @brief @ref Expression representing a scalar polynomial.
/// @ingroup fem
template <concepts::Numeric TValue, int PolynomialOrder = -1>
class PolynomialView : public ExpressionTraits<TValue> {
private:
    constexpr static inline int Order = PolynomialOrder;

    constexpr static inline bool hasStaticCoefficients = (0<= Order);

    constexpr static inline unsigned coefficientCount = hasStaticCoefficients
        ? static_cast<unsigned>(PolynomialOrder) + 1u
        : 0u;

public:
    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = std::conditional_t<
        hasStaticCoefficients,
        PolynomialView<TValue,std::max(0,PolynomialOrder-1)>,
        PolynomialView<TValue,-1>
    >;

    using Coefficients = std::conditional_t<
        hasStaticCoefficients,
        std::array<TValue,coefficientCount>,
        DynamicArray<TValue>>;

    constexpr PolynomialView() noexcept = default;

    PolynomialView(ConstSpan coefficients) noexcept
    requires (!hasStaticCoefficients);

    void evaluate(ConstSpan in, Span out) const
    requires (!hasStaticCoefficients);

    static constexpr unsigned size() noexcept
    requires (!hasStaticCoefficients);

    Derivative makeDerivative(Span buffer) const
    requires (!hasStaticCoefficients);

    ConstSpan coefficients() const noexcept
    requires (!hasStaticCoefficients);

    constexpr PolynomialView(std::span<const TValue,coefficientCount> coefficients) noexcept
    requires hasStaticCoefficients;

    constexpr void evaluate(ConstSpan in, Span out) const
    requires hasStaticCoefficients;

    static constexpr unsigned size() noexcept
    requires hasStaticCoefficients;

    constexpr Derivative makeDerivative(std::span<TValue,Derivative::coefficientCount> buffer) const noexcept
    requires hasStaticCoefficients;

    std::span<const TValue,coefficientCount> coefficients() const noexcept
    requires hasStaticCoefficients;

private:
    std::conditional_t<
        hasStaticCoefficients,
        std::span<const TValue,coefficientCount>,
        ConstSpan
    > _coefficients;
}; // class PolynomialView


/// @brief @ref Expression representing a scalar polynomial.
/// @ingroup fem
template <class TValue>
class Polynomial : public ExpressionTraits<TValue> {
public:
    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = Polynomial;

    using Coefficients = DynamicArray<TValue>;

    using View = PolynomialView<TValue>;

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
    Polynomial(ConstSpan coefficients);

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
struct io::GraphML::Serializer<maths::PolynomialView<TValue>> {
    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::PolynomialView<TValue>> rInstance);
}; // struct GraphML::Serializer<PolynomialView>


template <class TValue>
struct io::GraphML::Serializer<maths::Polynomial<TValue>> {
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
