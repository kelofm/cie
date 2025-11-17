#ifndef CIE_UTILS_STL_EXTENSION_HASH_HPP
#define CIE_UTILS_STL_EXTENSION_HASH_HPP

// --- STL Includes ---
#include <concepts> // default_initializable, equality_comparable
#include <functional> // std::equal_to


namespace cie::utils {


template <class T>
struct Hash
{Hash() = delete;};


template <class T>
struct Equal {Equal() = delete;};


template <class T>
requires std::default_initializable<std::hash<T>>
struct Hash<T>
{
    std::size_t operator()(const T& rInstance) const noexcept;
};


template <class TFirst, class TSecond>
requires (std::default_initializable<std::hash<TFirst>>
          && std::default_initializable<std::hash<TSecond>>)
struct Hash<std::pair<TFirst,TSecond>>
{
    std::size_t operator()(const std::pair<TFirst,TSecond>& rInstance) const noexcept;
}; // Hash<pair>


template <class TFirst, class TSecond>
requires (  !(std::default_initializable<std::hash<TFirst>>
              && std::default_initializable<std::hash<TSecond>>)
          && (std::default_initializable<Hash<TFirst>>
              && std::default_initializable<Hash<TSecond>>))
struct Hash<std::pair<TFirst,TSecond>>
{
    std::size_t operator()(const std::pair<TFirst,TSecond>& rInstance) const noexcept;
}; // Hash<pair>


template <class T>
requires std::default_initializable<std::equal_to<T>>
struct Equal<T>
{
    bool operator()(const T& rLeft, const T& rRight) const noexcept;
};


template <class TFirst, class TSecond>
requires (std::equality_comparable<TFirst>
          && std::equality_comparable<TSecond>)
struct Equal<std::pair<TFirst,TSecond>>
{
    bool operator()(const std::pair<TFirst,TSecond>& left,
                    const std::pair<TFirst,TSecond>& right) const noexcept;
}; // struct Equal


} // namespace cie::utils

#include "packages/stl_extension/impl/Hash_impl.hpp"

#endif
