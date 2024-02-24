#ifndef CIE_UTILS_STL_EXTENSION_HASH_HPP
#define CIE_UTILS_STL_EXTENSION_HASH_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"


namespace cie::utils {


template <class T>
struct Hash : private std::hash<T>
{
    std::size_t operator()(const T& r_value) const noexcept;
};


template <class T>
struct Equal {};


template <concepts::Integer TFirst, concepts::Integer TSecond>
struct Hash<std::pair<TFirst,TSecond>>
{
    std::size_t operator()(std::pair<TFirst,TSecond> value) const noexcept;
}; // struct Hash


template <concepts::Integer TFirst, concepts::Integer TSecond>
struct Equal<std::pair<TFirst,TSecond>>
{
    bool operator()(std::pair<TFirst,TSecond> left,
                    std::pair<TFirst,TSecond> right) const noexcept;
}; // struct Equaö


} // namespace cie::utils

#include "packages/stl_extension/impl/Hash_impl.hpp"

#endif
