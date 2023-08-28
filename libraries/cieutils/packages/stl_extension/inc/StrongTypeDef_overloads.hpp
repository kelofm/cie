#ifndef CIE_UTILS_STRONG_TYPEDEF_OVERLOADS_HPP
#define CIE_UTILS_STRONG_TYPEDEF_OVERLOADS_HPP

// --- Internal Includes ---
#include "packages/stl_extension/inc/StrongTypeDef.hpp"


namespace cie::utils {


/* --- LOGICAL OPERATORS --- */

template <class T, class Tag>
bool operator==(const StrongTypeDef<T,Tag>& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T, class Tag>
bool operator==(const StrongTypeDef<T,Tag>& r_lhs,
                const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator==(const T& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;


template <class T, class Tag>
bool operator!=(const StrongTypeDef<T,Tag>& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T, class Tag>
bool operator!=(const StrongTypeDef<T,Tag>& r_lhs,
                const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator!=(const T& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;


template <class T, class Tag>
bool operator<(const StrongTypeDef<T,Tag>& r_lhs,
               const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T,class Tag>
bool operator<(const StrongTypeDef<T,Tag>& r_lhs,
               const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator<(const T& r_lhs,
               const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T, class Tag>
bool operator>(const StrongTypeDef<T,Tag>& r_lhs,
               const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T,class Tag>
bool operator>(const StrongTypeDef<T,Tag>& r_lhs,
               const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator>(const T& r_lhs,
               const StrongTypeDef<T,Tag>& r_rhs) noexcept;


template <class T, class Tag>
bool operator<=(const StrongTypeDef<T,Tag>& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T,class Tag>
bool operator<=(const StrongTypeDef<T,Tag>& r_lhs,
                const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator<=(const T& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T, class Tag>
bool operator>=(const StrongTypeDef<T,Tag>& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;

template <class T,class Tag>
bool operator>=(const StrongTypeDef<T,Tag>& r_lhs,
                const T& r_rhs) noexcept;

template <class T, class Tag>
bool operator>=(const T& r_lhs,
                const StrongTypeDef<T,Tag>& r_rhs) noexcept;


/* --- ARITHMETIC OPERATORS --- */

/// @todo implement arithmetic operators for @ref StrongTypeDef


/* --- STREAM OPERATORS --- */

/// @todo implement stream operators for @ref StrongTypeDef


} // namespace cie::utils

#include "packages/stl_extension/impl/StrongTypeDef_overloads_impl.hpp"

#endif
