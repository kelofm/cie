#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/types/inc/tags.hpp"

// --- STL Includes ---
#include <iosfwd>


namespace cie::io {


/// @ingroup cieutils
struct Traits {
    using SerializerStream = std::ostream;

    using DeserializerStream = std::istream;
}; // struct Traits


struct TriviallySerializableBase {};


} // namespace cie::io


namespace cie {


/** @brief Concept for trivially serializable types.
 *  @details Trivial serialization means different thins depending
 *           on the type of serialization. For binary serialization,
 *           it means that the serialized output is identical to the
 *           unserialized binary representation. For text serialization,
 *           it means that the object can be directly serialized to an
 *           std::ostream and deserialized from an std::istream via
 *           operator<< and operator>> respectively.
 *  @note Trivial serialization implies trivial deserialization as well.
 *  @ingroup cieutils
 */
template <class T>
concept TriviallySerializable
=   std::integral<T>
    || std::floating_point<T>
    || std::is_base_of_v<cie::io::TriviallySerializableBase,T>;


/// @ingroup cieutils
template <class T>
concept BinarySerializable
= TriviallySerializable<T>
|| requires (Ref<const T> rInstance, Ref<io::Traits::SerializerStream> rStream) {
    {rInstance.serialize(rStream, tags::Binary())};
}; // concept BinarySerializable


/// @ingroup cieutils
template <class T>
concept TextSerializable
= TriviallySerializable<T>
|| requires (Ref<const T> rInstance, Ref<io::Traits::SerializerStream> rStream) {
    {rInstance.serialize(rStream, tags::Text())};
}; // concept TextSerializable


/// @ingroup cieutils
template <class T, class TTag = tags::Null>
concept Serializable
  = (std::is_same_v<TTag,tags::Null> && (BinarySerializable<T> || TextSerializable<T>))
    || (std::is_same_v<TTag,tags::Binary> && BinarySerializable<T>)
    || (std::is_same_v<TTag,tags::Text> && TextSerializable<T>);


/// @ingroup cieutils
template <class T>
concept NonTriviallySerializable
= !TriviallySerializable<T> && Serializable<T>;


/** @brief Concept for trivially deserializable types.
 *  @details Trivial deserialization means different things depending
 *           on the type of deserialization. For binary deserialization,
 *           it means that the deserialized output is identical to the
 *           serialized binary representation. For text serialization,
 *           it means that the object can be directly deserialized from an
 *           std::ostream and serialized to an std::istream via
 *           operator<< and operator>> respectively.
 *  @note Trivial deserialization implies trivial serialization as well.
 *  @ingroup cieutils
 */
template <class T>
concept TriviallyDeserializable
= TriviallySerializable<T>;


/// @ingroup cieutils
template <class T, class TAllocator>
concept BinaryDeserializable
= TriviallyDeserializable<T>
|| requires (
    Ref<io::Traits::DeserializerStream> rStream,
    Ref<T> rInstance,
    Ref<TAllocator> rAllocator) {
        {
            T::deserialize(
                rStream,
                rInstance,
                rAllocator,
                tags::Binary())
        };
};


/// @ingroup cieutils
template <class T, class TAllocator>
concept TextDeserializable
= TriviallyDeserializable<T>
|| requires (
    Ref<io::Traits::DeserializerStream> rStream,
    Ref<T> rInstance,
    Ref<TAllocator> rAllocator) {
        {
            T::deserialize(
                rStream,
                rInstance,
                rAllocator,
                tags::Text())
        };
};


/// @ingroup cieutils
template <class T, class TAllocator, class TTag = tags::Null>
concept Deserializable
  = (std::is_same_v<TTag,tags::Null> && (BinaryDeserializable<T,TAllocator> || TextDeserializable<T,TAllocator>))
    || (std::is_same_v<TTag,tags::Binary> && BinaryDeserializable<T,TAllocator>)
    || (std::is_same_v<TTag,tags::Text> && TextDeserializable<T,TAllocator>);


/// @ingroup cieutils
template <class T, class TAllocator>
concept NonTriviallyDeserializable
= !TriviallyDeserializable<T> && Deserializable<T,TAllocator>;


} // namespace cie
