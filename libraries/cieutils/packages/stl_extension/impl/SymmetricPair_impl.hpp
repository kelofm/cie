#ifndef CIE_UTILS_SYMMETRIC_PAIR_IMPL_HPP
#define CIE_UTILS_SYMMETRIC_PAIR_IMPL_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/SymmetricPair.hpp"


namespace cie::utils {


template <class TValue>
std::size_t Hash<SymmetricPair<TValue>>::operator()(Ref<const SymmetricPair<TValue>> r_pair) const noexcept
{
    const Hash<TValue> hash;
    return hash(r_pair.first) + hash(r_pair.second); // <== yes, this is a shitty idea @todo @fixme
}


template <class TValue>
bool operator==(Ref<const SymmetricPair<TValue>> r_lhs,
                Ref<const SymmetricPair<TValue>> r_rhs) noexcept
{
    return (
        (r_lhs.first == r_rhs.first)
        &&
        (r_lhs.second == r_rhs.second)
    ) || (
        (r_lhs.first == r_rhs.second)
        &&
        (r_lhs.second == r_rhs.first)
    );
}


template <class TValue>
bool operator!=(Ref<const SymmetricPair<TValue>> r_lhs,
                Ref<const SymmetricPair<TValue>> r_rhs) noexcept
{
    return (
        (r_lhs.first != r_rhs.first)
        ||
        (r_lhs.second != r_rhs.second)
    ) && (
        (r_lhs.first != r_rhs.second)
        ||
        (r_lhs.second != r_rhs.first)
    );
}


} // namespace cie::utils


#endif
