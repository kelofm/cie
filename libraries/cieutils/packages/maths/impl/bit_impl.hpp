#ifndef CIE_UTILS_BIT_OPERATIONS_IMPL_HPP
#define CIE_UTILS_BIT_OPERATIONS_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/bit.hpp"


namespace cie::utils {


template <concepts::Integer TInteger, concepts::Integer TIndex>
inline TInteger
flipBit(TInteger integer, TIndex bitIndex)
{
    return integer ^ (1 << bitIndex);
}


template <concepts::Integer TInteger, concepts::Integer TIndex>
bool getBit( TInteger integer, TIndex bitIndex)
{
    return (integer & (1 << bitIndex)) >> bitIndex;
}


template <concepts::UnsignedInteger TInteger>
std::uint8_t getMSBIndex(TInteger integer) noexcept
{
    std::uint8_t index = 0;
    while(integer >>= 1) {++index;}
    return index;
}


template <concepts::UnsignedInteger TInteger>
std::uint8_t getNumberOfTrailingZeros(TInteger integer) noexcept
{
    std::uint8_t count = 0;
    TInteger mask = 1;
    while (!(integer & mask) && mask) {mask <<= 1; ++count;}
    return count;
}


} // namespace cie::utils


#endif