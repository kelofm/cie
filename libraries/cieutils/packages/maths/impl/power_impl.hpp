#ifndef CIE_UTILS_MATHS_POWER_IMPL_HPP
#define CIE_UTILS_MATHS_POWER_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/assertions.hpp"


namespace cie {


template <concepts::Integer TBase, concepts::Integer TExponent>
constexpr TBase intPow(TBase base, TExponent exponent)
{
    // Negative exponents for integer bases yield non-integer results!
    assert(exponent >= 0);

    if (exponent == 0)
        return 1;
    else if (exponent == 1)
        return base;

    TBase output = intPow(base, exponent/2);

    if (exponent%2 == 0)
        return output * output;
    else
        return output * output * base;
}


} // namespace cie


#endif