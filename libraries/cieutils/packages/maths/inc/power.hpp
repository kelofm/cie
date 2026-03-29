#ifndef CIE_UTILS_MATHS_POWER_HPP
#define CIE_UTILS_MATHS_POWER_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie {


/// @ingroup cieutils
template <concepts::Integer TBase, concepts::Integer TExponent>
constexpr TBase intPow(TBase base, TExponent exponent);


} // namespace cie

#include "packages/maths/impl/power_impl.hpp"

#endif
