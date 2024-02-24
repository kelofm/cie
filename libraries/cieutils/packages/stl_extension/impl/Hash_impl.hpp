#ifndef CIE_UTILS_STL_EXTENSION_HASH_IMPL_HPP
#define CIE_UTILS_STL_EXTENSION_HASH_IMPL_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/Hash.hpp"

// --- STL Includes ---
#include <functional> // hash


namespace cie::utils {


template <class T>
std::size_t Hash<T>::operator()(const T& r_value) const noexcept
{
    return std::hash<T>().operator()(r_value);
}


template <concepts::Integer TFirst, concepts::Integer TSecond>
std::size_t Hash<std::pair<TFirst,TSecond>>::operator()(std::pair<TFirst,TSecond> value) const noexcept
{
    const std::size_t tmp = std::hash<TFirst>()(value.first);
    return tmp ^ (std::hash<cie::Size>()(value.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
}


template <concepts::Integer TFirst, concepts::Integer TSecond>
bool Equal<std::pair<TFirst,TSecond>>::operator()(std::pair<TFirst,TSecond> left,
                                                  std::pair<TFirst,TSecond> right) const noexcept
{
    return left.first == right.first && left.second == right.second;
}


} // namespace cie::utils


#endif
