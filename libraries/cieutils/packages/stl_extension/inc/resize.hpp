#ifndef CIE_CIEUTILS_CONCEPTS_UTILITY_HPP
#define CIE_CIEUTILS_CONCEPTS_UTILITY_HPP

// --- Internal Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <cassert>


namespace cie::utils {


/// @addtogroup cieutils
template <concepts::ResizableContainer ContainerType>
inline void resize(ContainerType& container, Size size)
{container.resize(size);}


/// @addtogroup cieutils
template <concepts::StaticContainer ContainerType>
inline void resize(ContainerType& container, Size size )
{CIE_CHECK(container.size() == size, "Attempt to resize static container")}


/// @addtogroup cieutils
template <concepts::ReservableContainer ContainerType>
inline void reserve(ContainerType& container, Size size)
{container.reserve(size);}


/// @addtogroup cieutils
template <concepts::StaticContainer ContainerType>
inline void reserve(ContainerType& container, Size size)
{CIE_CHECK(container.size() >= size, "Attempt to reserve static container")}


} // namespace cie::utils

#endif