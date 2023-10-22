#ifndef CIE_UTILS_SYMMETRIC_PAIR_HPP
#define CIE_UTILS_SYMMETRIC_PAIR_HPP

// --- Utiltiy Includes ---
#include "packages/stl_extension/inc/Hash.hpp"

// --- STL Includes ---
#include <functional> // hash


namespace cie::utils {


/// @brief Identical to std::pair but overloads comparison and hash operators to disregard member order.
template <class TValue>
struct SymmetricPair
{
    TValue first;

    TValue second;
};


template <class TValue>
struct Hash<SymmetricPair<TValue>>
{
    /// @todo The current implementation probably makes hash tables perform horribly => need a better implementation.
    std::size_t operator()(Ref<const SymmetricPair<TValue>> r_pair) const noexcept;
};


template <class TValue>
bool operator==(Ref<const SymmetricPair<TValue>> r_lhs,
                Ref<const SymmetricPair<TValue>> r_rhs) noexcept;


template <class TValue>
bool operator!=(Ref<const SymmetricPair<TValue>> r_lhs,
                Ref<const SymmetricPair<TValue>> r_rhs) noexcept;


} // namespace cie::utils


#include "packages/stl_extension/impl/SymmetricPair_impl.hpp"

#endif
