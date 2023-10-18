#ifndef CIE_UTILS_MATHS_BASE_HPP
#define CIE_UTILS_MATHS_BASE_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <array>


namespace cie::utils {


/// @addtogroup cieutils
template <Size NDigits, concepts::Integer NT>
constexpr void
baseN(NT base10,
      NT base,
      StaticArray<NT,NDigits>& output);


/// @addtogroup cieutils
template <Size NDigits, concepts::Integer NT>
constexpr StaticArray<NT,NDigits>
baseN(NT base10,
      NT base);


/// @addtogroup cieutils
template <Size NDigits, concepts::Integer NT>
constexpr NT
base10(const StaticArray<NT,NDigits>& base_N,
       NT base);


/// @addtogroup cieutils
template <Size NDigits, concepts::Integer NT>
constexpr NT
reinterpretBase(const StaticArray<NT,NDigits>& base_N,
                NT newBase);



} // namespace cie::utils

#include "packages/maths/impl/base_impl.hpp"

#endif
