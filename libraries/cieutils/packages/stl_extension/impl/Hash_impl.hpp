#ifndef CIE_UTILS_STL_EXTENSION_HASH_IMPL_HPP
#define CIE_UTILS_STL_EXTENSION_HASH_IMPL_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/Hash.hpp"

// --- STL Includes ---
#include <functional> // hash


namespace cie::utils {


template <class T>
requires std::default_initializable<std::hash<T>>
std::size_t Hash<T>::operator()(const T& rInstance) const noexcept
{
    return std::hash<T>()(rInstance);
}


template <class TFirst, class TSecond>
requires (std::default_initializable<std::hash<TFirst>>
          && std::default_initializable<std::hash<TSecond>>)
std::size_t Hash<std::pair<TFirst,TSecond>>::operator()(const std::pair<TFirst,TSecond>& rInstance) const noexcept
{
    const std::size_t tmp = std::hash<TFirst>()(rInstance.first);
    return tmp ^ (std::hash<TSecond>()(rInstance.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
}


template <class TFirst, class TSecond>
requires (std::default_initializable<Hash<TFirst>>
          && std::default_initializable<Hash<TSecond>>)
std::size_t Hash<std::pair<TFirst,TSecond>>::operator()(const std::pair<TFirst,TSecond>& rInstance) const noexcept
{
    const std::size_t tmp = Hash<TFirst>()(rInstance.first);
    return tmp ^ (Hash<TSecond>()(rInstance.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
}


template <class T>
requires std::default_initializable<std::equal_to<T>>
bool Equal<T>::operator()(const T& rLeft, const T& rRight) const noexcept
{
    return std::equal_to<T>()(rLeft, rRight);
}


template <class TFirst, class TSecond>
requires (std::equality_comparable<TFirst>
          && std::equality_comparable<TSecond>)
bool Equal<std::pair<TFirst,TSecond>>::operator()(const std::pair<TFirst,TSecond>& left,
                                                  const std::pair<TFirst,TSecond>& right) const noexcept
{
    return left.first == right.first && left.second == right.second;
}


} // namespace cie::utils


#endif
