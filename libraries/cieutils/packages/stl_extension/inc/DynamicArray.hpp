#ifndef CIE_UTILS_STL_EXTENSION_DYNAMIC_ARRAY_HPP
#define CIE_UTILS_STL_EXTENSION_DYNAMIC_ARRAY_HPP

// --- Internal Includes ---
#include "packages/stl_extension/inc/STLContainerBase.hpp"

// --- STL Includes ---
#include <vector>


namespace cie {


/// @ingroup cieutils
template <class TValue, class ...TArgs>
using DynamicArray = STLContainerBase<std::vector<TValue,TArgs...>>;


} // namespace cie


#endif
