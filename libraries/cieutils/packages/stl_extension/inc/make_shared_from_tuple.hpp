#ifndef CIE_UTILS_MAKE_SHARED_FROM_TUPLE_HPP
#define CIE_UTILS_MAKE_SHARED_FROM_TUPLE_HPP

// --- STL Includes ---
#include <tuple>
#include <memory>


namespace cie::utils {


/// @ingroup cieutils
template <class T, class ...Args>
std::shared_ptr<T> make_shared_from_tuple( std::tuple<Args...>& args );


} // namespace cie::utils

#include "packages/stl_extension/impl/make_shared_from_tuple_impl.hpp"

#endif
