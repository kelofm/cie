#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/StrongTypeDef.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"


namespace cie::fem {


/// @brief Vertex identifier type.
/// @class VertexID
/// @see StrongTypeDef
CIE_STRONG_TYPEDEF(unsigned, VertexID);

/// @brief Edge identifier type.
/// @class EdgeID
/// @see StrongTypeDef
CIE_STRONG_TYPEDEF(unsigned, EdgeID);


template <class T>
concept VertexLike
=  std::is_same_v<typename T::ID,VertexID>
&& requires (T& instance, const T& constInstance) {
    typename T::Data;
    {T(VertexID())}             -> std::same_as<T>;
    {constInstance.id()}        -> std::same_as<VertexID>;
    {instance.data()}           -> std::same_as<typename VoidSafe<typename T::Data>::Ref>;
    {constInstance.data()}      -> std::same_as<typename VoidSafe<const typename T::Data>::Ref>;
    {constInstance.edges()};
}; // concept VertexLike


template <class T>
concept EdgeLike
=  std::is_same_v<typename T::ID,EdgeID>
&& requires (T& instance, const T& constInstance) {
    typename T::Data;
    {T(EdgeID(), std::pair<VertexID,VertexID>())}   -> std::same_as<T>;
    {constInstance.id()}                            -> std::same_as<EdgeID>;
    {instance.data()}                               -> std::same_as<typename VoidSafe<typename T::Data>::Ref>;
    {constInstance.data()}                          -> std::same_as<typename VoidSafe<const typename T::Data>::Ref>;
    {constInstance.source()}                        -> std::same_as<VertexID>;
    {constInstance.target()}                        -> std::same_as<VertexID>;
    {constInstance.vertices()};
}; // concept VertexLike


template <class T>
concept GraphLike
=  VertexLike<typename T::Vertex>
&& EdgeLike<typename T::Edge>
&& requires (T& instance, const T& constInstance) {
    {instance.data()}                   -> std::same_as<typename VoidSafe<typename T::Data>::Ref>;
    {constInstance.data()}              -> std::same_as<typename VoidSafe<const typename T::Data>::Ref>;
    {instance.find(VertexID())}         -> std::same_as<OptionalRef<typename T::Vertex>>;
    {constInstance.find(VertexID())}    -> std::same_as<OptionalRef<const typename T::Vertex>>;
    {instance.find(EdgeID())}           -> std::same_as<OptionalRef<typename T::Edge>>;
    {constInstance.find(EdgeID())}      -> std::same_as<OptionalRef<const typename T::Edge>>;
}; // concept GraphLike


} // namespace cie::fem
