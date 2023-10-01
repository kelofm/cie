#ifndef CIE_UTILS_IO_TRAITS_HPP
#define CIE_UTILS_IO_TRAITS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/types/inc/tags.hpp"

// --- STL Includes ---
#include <iosfwd>


namespace cie::io {


struct Traits
{
    using SerializerStream = std::ostream;

    using DeserializerStream = std::istream;
}; // struct Traits


struct TriviallySerializable {
    inline constexpr static bool isTriviallySerializable = true;
};


} // namespace cie::io


namespace cie::concepts {


/** @brief Concept for trivially serializable types.
 *  @details Trivial serialization means different thins depending
 *           on the type of serialization. For binary serialization,
 *           it means that the serialized output is identical to the
 *           unserialized binary representation. For text serialization,
 *           it means that the object can be directly serialized to an
 *           std::ostream and deserialized from an std::istream via
 *           operator<< and operator>> respectively.
 *  @note Trivial serialization implies trivial deserialization as well.
 */
template <class T>
concept TriviallySerializable
=   std::integral<T>
    || std::floating_point<T>
    || T::isTriviallySerializable;


template <class T>
concept BinarySerializable
= TriviallySerializable<T>
|| requires (Ref<const T> r_instance, Ref<io::Traits::SerializerStream> r_stream)
{
    {r_instance.serialize(r_stream, tags::Binary())};
}; // concept BinarySerializable


template <class T>
concept TextSerializable
= TriviallySerializable<T>
|| requires (Ref<const T> r_instance, Ref<io::Traits::SerializerStream> r_stream)
{
    {r_instance.serialize(r_stream, tags::Text())};
}; // concept TextSerializable


template <class T, class TTag = tags::Null>
concept Serializable
  = (std::is_same_v<TTag,tags::Null> && (BinarySerializable<T> || TextSerializable<T>))
    || (std::is_same_v<TTag,tags::Binary> && BinarySerializable<T>)
    || (std::is_same_v<TTag,tags::Text> && TextSerializable<T>);


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
 */
template <class T>
concept TriviallyDeserializable
= TriviallySerializable<T>;


template <class T>
concept BinaryDeserializable
= TriviallyDeserializable<T>
|| requires (Ref<T> r_instance, Ref<io::Traits::DeserializerStream> r_stream)
{
    {r_instance.deserialize(r_stream, r_instance, tags::Binary())};
};


template <class T>
concept TextDeserializable
= TriviallyDeserializable<T>
|| requires (Ref<T> r_instance, Ref<io::Traits::DeserializerStream> r_stream)
{
    {r_instance.deserialize(r_stream, r_instance, tags::Text())};
};


template <class T, class TTag = tags::Null>
concept Deserializable
  = (std::is_same_v<TTag,tags::Null> && (BinaryDeserializable<T> || TextDeserializable<T>))
    || (std::is_same_v<TTag,tags::Binary> && BinaryDeserializable<T>)
    || (std::is_same_v<TTag,tags::Text> && TextDeserializable<T>);


template <class T>
concept NonTriviallyDeserializable
= !TriviallyDeserializable<T> && Deserializable<T>;


} // namespace cie::concepts


#endif
