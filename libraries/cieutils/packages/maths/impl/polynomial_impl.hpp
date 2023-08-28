#ifndef CIE_UTILS_MATHS_POLYNOMIAL_IMPL_HPP
#define CIE_UTILS_MATHS_POLYNOMIAL_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/polynomial.hpp"

// --- STL Includes ---
#include <cmath>
#include <limits>


namespace cie::utils {


template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
inline TValue evaluatePolynomialNaive(const TValue argument, TItBegin it_begin, const TItEnd it_end) noexcept
{
    TValue result = static_cast<TValue>(0);
    TValue power = static_cast<TValue>(1);
    for (; it_begin!=it_end; ++it_begin)
    {
        result += power * (*it_begin);
        power *= argument;
    }
    return result;
}


template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
inline TValue evaluatePolynomialHorner(const TValue argument, const TItBegin it_begin, TItEnd it_end) noexcept
{
    TValue result = static_cast<TValue>(0);

    // Reverse loop
    if (it_begin != it_end) [[likely]]
    {
        --it_end;
        do
        {
            result *= argument;
            result += *it_end;
        } while (it_end-- != it_begin);
    }

    return result;
}


namespace detail {
template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
inline TValue evaluatePolynomialHornerStabilizedNoCheck(TValue argument, TValue argumentInverse, TItBegin it_begin, const TItEnd it_end) noexcept
{
    TValue result = static_cast<TValue>(0);

    for (; it_begin!=it_end; ++it_begin)
    {
        result *= argumentInverse;
        result += *it_begin;
    }

    const unsigned polynomialOrder = std::distance(it_begin, it_end);
    if (polynomialOrder) [[likely]]
        result *= std::pow(argument, polynomialOrder - 1);

    return result;
}
} // namespace detail


template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
inline TValue evaluatePolynomialHornerStabilized(TValue argument, TItBegin it_begin, const TItEnd it_end) noexcept
{
    TValue result;
    constexpr const TValue threshold = std::numeric_limits<TValue>::is_integer ? static_cast<TValue>(1) : std::numeric_limits<TValue>::epsilon();
    if (threshold < std::abs(argument)) [[likely]]
    {
        const TValue argumentInverse = static_cast<TValue>(1) / argument;
        result = detail::evaluatePolynomialHornerStabilizedNoCheck(argument, argumentInverse, it_begin, it_end);
    }
    else
        result = static_cast<TValue>(0);

    return result;
}


template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
inline TValue evaluatePolynomialHornerCompound(TValue argument, TItBegin it_begin, TItEnd it_end) noexcept
{
    TValue result;
    if (static_cast<TValue>(1) < std::abs(argument))
    {
        const TValue argumentInverse = static_cast<TValue>(1) / argument;
        result = detail::evaluatePolynomialHornerStabilizedNoCheck(argument, argumentInverse, it_begin, it_end);
    }
    else
        result = evaluatePolynomialHorner(argument, it_begin, it_end);

    return result;
}


} // namespace cie::utils


#endif
