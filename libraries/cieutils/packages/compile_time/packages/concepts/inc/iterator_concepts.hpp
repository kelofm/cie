#ifndef CIE_CIEUTILS_ITERATOR_CONCEPTS_HPP
#define CIE_CIEUTILS_ITERATOR_CONCEPTS_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Select.hpp"

// --- STL Includes ---
#include <iterator>
#include <memory>


namespace cie::concepts {


/// @brief Iterator concept [optionally pointing to a specific type].
template <class ...TArgs>
concept Iterator
= !std::is_same_v<typename std::iterator_traits<typename ct::Select<TArgs...>::Front>::value_type,void>
  && (ct::PackSize<TArgs...> == 1 ? true : std::is_same_v<typename std::iterator_traits<typename ct::Select<TArgs...>::Front>::value_type,typename ct::Select<TArgs...>::Back>)
  && ct::PackSize<TArgs...> <= 2;


/// @brief Iterator over types that can be converted to @a TValue.
/// @tparam TIt: iterator type.
/// @tparam TValue: value type must be convertible to this type.
template <class TIt, class TValue>
concept WeakIterator
= Iterator<TIt> && std::convertible_to<typename std::iterator_traits<TIt>::value_type,TValue>;


} // namespace cie::concepts

#endif
