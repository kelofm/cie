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


template <concepts::Numeric TValue, int PolynomialOrder>
class Polynomial;


/// @brief @ref Expression representing a scalar polynomial.
/// @ingroup fem
template <concepts::Numeric TValue, int PolynomialOrder = -1>
class PolynomialView : public ExpressionTraits<TValue> {
private:
    constexpr static inline int Order = PolynomialOrder;

    constexpr static inline bool hasStaticCoefficients = (0<= Order);

public:
    constexpr static inline unsigned coefficientCount = hasStaticCoefficients
        ? static_cast<unsigned>(PolynomialOrder) + 1u
        : 0u;

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

    constexpr PolynomialView() noexcept;

    PolynomialView(ConstSpan coefficients) noexcept
    requires (!hasStaticCoefficients);

    void evaluate(ConstSpan in, Span out) const
    requires (!hasStaticCoefficients);

    static constexpr unsigned size() noexcept;

    Derivative makeDerivative(Span buffer) const
    requires (!hasStaticCoefficients);

    ConstSpan coefficients() const noexcept
    requires (!hasStaticCoefficients);

    constexpr PolynomialView(std::span<const TValue,coefficientCount> coefficients) noexcept
    requires (hasStaticCoefficients);

    constexpr void evaluate(ConstSpan in, Span out) const
    requires (hasStaticCoefficients);

    constexpr Derivative makeDerivative(std::span<TValue,Derivative::coefficientCount> buffer) const noexcept
    requires (hasStaticCoefficients);

    constexpr std::span<const TValue,coefficientCount> coefficients() const noexcept
    requires (hasStaticCoefficients);

private:
    friend class Polynomial<TValue,PolynomialOrder>;

    friend class Polynomial<TValue,PolynomialOrder+1>;

    std::conditional_t<
        hasStaticCoefficients,
        std::span<const TValue,coefficientCount>,
        ConstSpan
    > _coefficients;
}; // class PolynomialView


/// @brief @ref Expression representing a scalar polynomial.
/// @ingroup fem
template <concepts::Numeric TValue, int PolynomialOrder = -1>
class Polynomial : public ExpressionTraits<TValue> {
public:
    using View = PolynomialView<TValue,PolynomialOrder>;

    static inline constexpr bool hasStaticCoefficients = View::hasStaticCoefficients;

    static inline constexpr unsigned coefficientCount = View::coefficientCount;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = std::conditional_t<
        hasStaticCoefficients,
        Polynomial<TValue,std::max(0,PolynomialOrder-1)>,
        Polynomial<TValue,-1>>;

    using Coefficients = std::conditional_t<
        hasStaticCoefficients,
        std::array<TValue,coefficientCount>,
        DynamicArray<TValue>>;

public:
    /// @brief Uninitialized by default.
    constexpr Polynomial() noexcept = default;

    constexpr Polynomial(Polynomial&&) noexcept;

    Polynomial(Polynomial&& rRight)
    requires (!hasStaticCoefficients);

    constexpr Polynomial(Polynomial&& rRight)
    requires (hasStaticCoefficients);

    Polynomial(const Polynomial& rRight)
    requires (!hasStaticCoefficients);

    constexpr Polynomial(const Polynomial& rRight)
    requires (hasStaticCoefficients);

    /// @brief Construct from a container of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    constexpr Polynomial(RightRef<Coefficients> rCoefficients) noexcept;

    /// @brief Construct from a range of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    Polynomial(ConstSpan coefficients)
    requires (!hasStaticCoefficients);

    /// @brief Construct from a range of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    constexpr Polynomial(std::span<const TValue,coefficientCount> coefficients)
    requires (hasStaticCoefficients);

    Polynomial& operator=(Polynomial&&) noexcept
    requires (!hasStaticCoefficients);

    constexpr Polynomial& operator=(Polynomial&&) noexcept
    requires (hasStaticCoefficients);

    Polynomial& operator=(const Polynomial& rRight)
    requires (!hasStaticCoefficients);

    constexpr Polynomial& operator=(const Polynomial& rRight) noexcept
    requires (hasStaticCoefficients);

    void evaluate(ConstSpan in, Span out) const
    requires (!hasStaticCoefficients);

    constexpr void evaluate(ConstSpan in, Span out) const noexcept
    requires (hasStaticCoefficients);

    /// @brief Get the number of scalar components returned by @ref evaluate.
    constexpr static unsigned size() noexcept;

    /// @brief Construct the derivative of the @ref Polynomial.
    /// @details The returned object is also a @ref Polynomial,
    ///          irrespective of the current polynomial order.
    Derivative makeDerivative() const
    requires (!hasStaticCoefficients);

    /// @brief Construct the derivative of the @ref Polynomial.
    /// @details The returned object is also a @ref Polynomial,
    ///          irrespective of the current polynomial order.
    constexpr Derivative makeDerivative() const noexcept
    requires (hasStaticCoefficients);

    std::span<const TValue> coefficients() const noexcept
    requires (!hasStaticCoefficients);

    constexpr std::span<const TValue,coefficientCount> coefficients() const noexcept
    requires (hasStaticCoefficients);

private:
    template <concepts::Numeric T, int O>
    friend class Polynomial;

    Coefficients _coefficients;

    PolynomialView<TValue,PolynomialOrder> _wrapped;
}; // class Polynomial


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TValue, int PolynomialOrder>
struct io::GraphML::Serializer<maths::PolynomialView<TValue,PolynomialOrder>> {
    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::PolynomialView<TValue,PolynomialOrder>> rInstance);
}; // struct GraphML::Serializer<PolynomialView>


template <class TValue, int PolynomialOrder>
struct io::GraphML::Serializer<maths::Polynomial<TValue,PolynomialOrder>> {
    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::Polynomial<TValue,PolynomialOrder>> rInstance);
}; // struct GraphML::Serializer<Polynomial>


template <class TValue, int PolynomialOrder>
struct io::GraphML::Deserializer<maths::Polynomial<TValue,PolynomialOrder>>
    : public io::GraphML::DeserializerBase<maths::Polynomial<TValue,PolynomialOrder>>
{
    using io::GraphML::DeserializerBase<maths::Polynomial<TValue,PolynomialOrder>>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);

private:
    typename maths::Polynomial<TValue,PolynomialOrder>::Coefficients _coefficients;
}; // struct GraphML::Deserializer<Polynomial>


} // namespace cie::fem::io


#include "packages/maths/impl/Polynomial_impl.hpp"
