#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- STL Includes ---
#include <span>


namespace cie::utils {


enum class PolynomialEvaluation {
    Naive,
    Horner,
    HornerStabilized,
    HornerCompound
}; // enum class PolynomialEvaluation


/// @brief Evaluate a polynomial.
/// @tparam TEvaluation Evaluation strategy.
/// @tparam TValue Value type of the argument and the result.
/// @tparam TCoefficient Value type of the coefficients.
/// @param argument Position at which to evaluate the polynomial.
/// @param coefficients Array of coefficients in ascending monomial order.
/// @note The coefficients are expected to be sorted with respect to their order in the polynomial
///       (increasing order, starting with order 0 and including null value coefficients).
/// @warning This function does not check for overflows.
/// @ingroup cieutils
template <
    PolynomialEvaluation TEvaluation = PolynomialEvaluation::Horner,
    concepts::Numeric TValue,
    concepts::Numeric TCoefficient,
    std::size_t CoefficientCount>
requires (CoefficientCount == std::dynamic_extent)
[[nodiscard]] TValue evaluatePolynomial(
    TValue argument,
    std::span<const TCoefficient,CoefficientCount> coefficients) noexcept;


/// @brief Evaluate a polynomial.
/// @tparam TEvaluation Evaluation strategy.
/// @tparam TValue Value type of the argument and the result.
/// @tparam TCoefficient Value type of the coefficients.
/// @param argument Position at which to evaluate the polynomial.
/// @param coefficients Array of coefficients in ascending monomial order.
/// @note The coefficients are expected to be sorted with respect to their order in the polynomial
///       (increasing order, starting with order 0 and including null value coefficients).
/// @warning This function does not check for overflows.
/// @ingroup cieutils
template <
    PolynomialEvaluation TEvaluation = PolynomialEvaluation::Horner,
    concepts::Numeric TValue,
    concepts::Numeric TCoefficient,
    std::size_t CoefficientCount>
requires (CoefficientCount != std::dynamic_extent)
[[nodiscard]] constexpr TValue evaluatePolynomial(
    TValue argument,
    std::span<const TCoefficient,CoefficientCount> coefficients) noexcept;


} // namespace cie::utils

#include "packages/maths/impl/polynomial_evaluation_impl.hpp"
