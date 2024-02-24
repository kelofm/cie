#ifndef CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SIZE_HPP
#define CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SIZE_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"



namespace cie::ct {


/// @addtogroup cieutils
template <class ...Ts>
inline constexpr Size PackSize = sizeof...(Ts);


/// @addtogroup cieutils
template <class ...Ts>
inline constexpr bool Empty = PackSize<Ts...> == 0;


} // namespace cie::ct


#endif