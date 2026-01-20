#pragma once

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"

// --- STL Includes ---
#include <cmath>
#include <limits>


namespace cie::utils {


template <
    PolynomialEvaluation TEvaluation,
    concepts::Numeric TValue,
    concepts::Numeric TCoefficient,
    std::size_t CoefficientCount>
requires (CoefficientCount == std::dynamic_extent)
[[nodiscard]] TValue evaluatePolynomial(TValue argument, std::span<const TCoefficient,CoefficientCount> coefficients) noexcept {
    const TCoefficient* it = coefficients.data();
    const TCoefficient* pEnd = it + coefficients.size();

    if constexpr (TEvaluation == PolynomialEvaluation::Naive) {
        TValue result = static_cast<TValue>(0);
        TValue power = static_cast<TValue>(1);
        for (; it!=pEnd; ++it) {
            result += power * (*it);
            power *= argument;
        }
        return result;
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::Horner) {
        TValue result = static_cast<TValue>(0);
        for (auto it=coefficients.rbegin(); it!=coefficients.rend(); ++it) {
            result *= argument;
            result += *it;
        }
        return result;
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::HornerStabilized) {
        constexpr const TValue threshold = std::numeric_limits<TValue>::is_integer
            ? static_cast<TValue>(1)
            : std::numeric_limits<TValue>::epsilon();
        if (threshold < std::abs(argument)) [[likely]] {
            const TValue argumentInverse = static_cast<TValue>(1) / argument;
            TValue result = static_cast<TValue>(0);

            for (; it!=pEnd; ++it) {
                result *= argumentInverse;
                result += *it;
            }

            const unsigned polynomialOrder = std::distance(it, pEnd);
            if (polynomialOrder) [[likely]]
                result *= std::pow(argument, polynomialOrder - 1);

            return result;
        } else {
            return 0;
        }
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::HornerCompound) {
        if (static_cast<TValue>(1) < std::abs(argument)) {
            const TValue argumentInverse = static_cast<TValue>(1) / argument;
            TValue result = static_cast<TValue>(0);

            for (; it!=pEnd; ++it) {
                result *= argumentInverse;
                result += *it;
            }

            const unsigned polynomialOrder = std::distance(it, pEnd);
            if (polynomialOrder) [[likely]]
                result *= std::pow(argument, polynomialOrder - 1);

            return result;
        } else {
            TValue result = static_cast<TValue>(0);
            for (auto it=coefficients.rbegin(); it!=coefficients.rend(); ++it) {
                result *= argument;
                result += *it;
            }
            return result;
        }
    }
}


template <
    PolynomialEvaluation TEvaluation,
    concepts::Numeric TValue,
    concepts::Numeric TCoefficient,
    std::size_t CoefficientCount>
requires (CoefficientCount != std::dynamic_extent)
[[nodiscard]] constexpr TValue evaluatePolynomial(TValue argument, std::span<const TCoefficient,CoefficientCount> coefficients) noexcept {
    [[maybe_unused]] const TCoefficient* it = coefficients.data();
    [[maybe_unused]] const TCoefficient* pEnd = it + coefficients.size();

    if constexpr (TEvaluation == PolynomialEvaluation::Naive) {
        TValue result = static_cast<TValue>(0);
        TValue power = static_cast<TValue>(1);
        for (; it!=pEnd; ++it) {
            result += power * (*it);
            power *= argument;
        }
        return result;
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::Horner) {
        TValue result = static_cast<TValue>(0);
        for (auto it=coefficients.rbegin(); it!=coefficients.rend(); ++it) {
            result *= argument;
            result += *it;
        }
        return result;
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::HornerStabilized) {
        constexpr const TValue threshold = std::numeric_limits<TValue>::is_integer
            ? static_cast<TValue>(1)
            : std::numeric_limits<TValue>::epsilon();
        if (threshold < std::abs(argument)) [[likely]] {
            const TValue argumentInverse = static_cast<TValue>(1) / argument;
            TValue result = static_cast<TValue>(0);

            for (; it!=pEnd; ++it) {
                result *= argumentInverse;
                result += *it;
            }

            const unsigned polynomialOrder = std::distance(it, pEnd);
            if (polynomialOrder) [[likely]]
                result *= std::pow(argument, polynomialOrder - 1);

            return result;
        } else {
            return 0;
        }
    }


    else if constexpr (TEvaluation == PolynomialEvaluation::HornerCompound) {
        if (static_cast<TValue>(1) < std::abs(argument)) {
            const TValue argumentInverse = static_cast<TValue>(1) / argument;
            TValue result = static_cast<TValue>(0);

            for (; it!=pEnd; ++it) {
                result *= argumentInverse;
                result += *it;
            }

            const unsigned polynomialOrder = std::distance(it, pEnd);
            if (polynomialOrder) [[likely]]
                result *= std::pow(argument, polynomialOrder - 1);

            return result;
        } else {
            TValue result = static_cast<TValue>(0);
            for (auto it=coefficients.rbegin(); it!=coefficients.rend(); ++it) {
                result *= argument;
                result += *it;
            }
            return result;
        }
    }
}


} // namespace cie::utils
