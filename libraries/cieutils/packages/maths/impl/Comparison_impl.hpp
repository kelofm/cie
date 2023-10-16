#ifndef CIE_UTILS_COMPARISON_IMPL_HPP
#define CIE_UTILS_COMPARISON_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"

// --- STL Includes ---
#include <limits>
#include <cmath>


namespace cie::utils {


template <class TFloat>
FloatComparison<TFloat>::FloatComparison() noexcept
    : FloatComparison(std::numeric_limits<TFloat>::min(),
                      std::numeric_limits<TFloat>::epsilon())
{
}


template <class TFloat>
FloatComparison<TFloat>::FloatComparison(const TFloat absoluteTolerance,
                                         const TFloat relativeTolerance) noexcept
    : _absoluteTolerance(absoluteTolerance),
      _relativeTolerance(relativeTolerance)
{
}


template <class TFloat>
bool FloatComparison<TFloat>::equal(const TFloat left,
                                    const TFloat right) const noexcept
{
    const TFloat norm = std::min(std::abs(left) + std::abs(right),
                                 std::numeric_limits<TFloat>::max());
    return std::abs(left - right) < std::max(_absoluteTolerance, _relativeTolerance * norm);
}


} // namespace cie::utils


#endif
