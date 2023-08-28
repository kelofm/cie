#ifndef CIE_UTILS_STL_EXTENSION_GET_REFERENCE_HPP
#define CIE_UTILS_STL_EXTENSION_GET_REFERENCE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::utils {


template <class T>
requires concepts::Pointer<T>
inline typename std::pointer_traits<typename std::decay<T>::type>::element_type&
getRef(T& r_p)
{
    CIE_CHECK_POINTER(r_p)
    return *r_p;
}


template <class T>
requires concepts::Pointer<T>
inline const typename std::pointer_traits<typename std::decay<T>::type>::element_type&
getRef(const T& r_p)
{
    CIE_CHECK_POINTER(r_p)
    return *r_p;
}


template <class T>
requires concepts::NonPointer<T>
inline T&
getRef(T& r)
{
    return r;
}


template <class T>
requires concepts::NonPointer<T>
inline const T&
getRef(const T& r)
{
    return r;
}


} // namespace cie::utils


#endif