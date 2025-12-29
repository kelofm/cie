#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::fem::concepts {


template <class T>
concept Element =
   std::is_same_v<std::remove_cvref_t<decltype(T::Dimension)>,unsigned>
&& ::cie::concepts::Numeric<typename T::Value>
&& requires (const T& rConstInstance,
             std::span<const typename Kernel<T::Dimension,typename T::Value>::LocalCoordinate,T::Dimension> constLocalSpan,
             std::span<typename Kernel<T::Dimension,typename T::Value>::LocalCoordinate,T::Dimension> localSpan,
             std::span<const typename Kernel<T::Dimension,typename T::Value>::GlobalCoordinate,T::Dimension> constGlobalSpan,
             std::span<typename Kernel<T::Dimension,typename T::Value>::GlobalCoordinate,T::Dimension> globalSpan) {
    {rConstInstance.transform(constLocalSpan, globalSpan)}  -> std::same_as<void>;      // <== transform from local to global space
    {rConstInstance.transform(constGlobalSpan, localSpan)}  -> std::same_as<void>;      // <== transform from global to local space
    {rConstInstance.ID()}                                   -> std::same_as<VertexID>;
    {rConstInstance.ansatzSpaceID()}                        -> ::cie::concepts::UnsignedInteger;
}; // concept Element


template <class T>
concept Cell = Element<T>;


} // namespace cie::fem::concepts
