#ifndef CIE_UTILS_MATHS_POWER_HPP
#define CIE_UTILS_MATHS_POWER_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie {


/// @addtogroup cieutils
template <concepts::Integer BaseType, concepts::Integer ExponentType>
constexpr BaseType intPow(BaseType base, ExponentType exponent);


} // namespace cie

#include "packages/maths/impl/power_impl.hpp"

#endif
