#ifndef CIE_UTILS_BIT_OPERATIONS_HPP
#define CIE_UTILS_BIT_OPERATIONS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <cstdint>


namespace cie::utils {


/// @ingroup cieutils
template <concepts::Integer TInteger, concepts::Integer TIndex = Size>
TInteger flipBit(TInteger integer, TIndex bitIndex);


/// @ingroup cieutils
template <concepts::Integer TInteger, concepts::Integer TIndex = Size>
bool getBit(TInteger integer, TIndex bitIndex);


/// @brief Get the index of the most significat set bit (counting from LSB)
/// @note Returns 0 if the input is 0.
/// @ingroup cieutils
template <concepts::UnsignedInteger TInteger>
constexpr std::uint8_t getMSBIndex(TInteger integer) noexcept;


/// @brief Find the number of trailing zeros of the input integer's binary representation.
/// @todo Add compiler-specifics.
/// @ingroup cieutils
template <concepts::UnsignedInteger TInteger>
std::uint8_t getNumberOfTrailingZeros(TInteger integer) noexcept;


} // namespace cie::utils

#include "packages/maths/impl/bit_impl.hpp"

#endif