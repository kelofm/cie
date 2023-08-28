#ifndef CIE_CONCEPTS_LOGICAL_HPP
#define CIE_CONCEPTS_LOGICAL_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Select.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Range.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"


namespace cie::concepts {


/// True if the last parameter matches any of the previous ones.
template <class T, class ...Ts>
concept Any
= ct::Match<T>::template Any<Ts...>;


/// True if the last parameter doesn't match any of the previous ones.
template <class T, class ...Ts>
concept None
= !Any<T,Ts...>;


} // namespace cie::concepts


#endif