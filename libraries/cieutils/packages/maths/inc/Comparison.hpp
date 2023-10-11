#ifndef CIE_UTILS_COMPARISON_HPP
#define CIE_UTILS_COMPARISON_HPP

// --- STL Includes ---
#include <concepts>


namespace cie::utils {


template <class TFloat>
class FloatComparison
{
public:
    FloatComparison() noexcept;

    FloatComparison(const TFloat absoluteTolerance,
                    const TFloat relativeTolerance) noexcept;

    bool equal(const TFloat left,
               const TFloat right) const noexcept;

private:
    TFloat _absoluteTolerance;

    TFloat _relativeTolerance;
}; // class Comparison


template <class T>
class Comparison {};


template <std::floating_point TFloat>
class Comparison<TFloat> : public FloatComparison<TFloat>
{
public:
    using FloatComparison<TFloat>::FloatComparison;
};


} // namespace cie::utils

#include "packages/maths/impl/Comparison_impl.hpp"

#endif
