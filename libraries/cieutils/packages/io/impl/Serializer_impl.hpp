#ifndef CIE_UTILS_SERIALIZER_IMPL_HPP
#define CIE_UTILS_SERIALIZER_IMPL_HPP

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"
#include "packages/io/inc/Traits.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::io {


template <>
template <concepts::TriviallySerializable T>
inline void Serializer<tags::Binary>::serialize(Ref<SerializerStream> r_stream, T object)
{
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state before trivial serialization")
    Serializer::serializeImpl(r_stream, object);
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state after trivial serialization")
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <concepts::NonTriviallySerializable T>
inline void Serializer<TTag>::serialize(Ref<SerializerStream> r_stream, Ref<const T> r_object)
{
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state before serialization")
    Serializer::serializeImpl(r_stream, r_object);
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state after serialization")
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <class T>
requires concepts::Serializable<T,TTag>
inline void Serializer<TTag>::serialize(Ref<SerializerStream> r_stream,
                                        Ptr<const T> begin,
                                        Size numberOfItems)
{
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state before array serialization")
    const Ptr<const T> end = begin + numberOfItems;
    for (; begin<end; ++begin)
        Serializer::serialize(r_stream, *begin);
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state after array serialization")
}


template <>
template <concepts::Deserializable T>
inline void Serializer<tags::Binary>::deserialize(Ref<DeserializerStream> r_stream, Ref<T> r_output)
{
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state before trivial deserialization")
    Serializer::deserializeImpl(r_stream, r_output);
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state after trivial deserialization")
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <concepts::Deserializable<TTag> T>
inline void Serializer<TTag>::deserialize(Ref<DeserializerStream> r_stream,
                                          Ptr<T> begin,
                                          Size numberOfItems)
{
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state before array deserialization")
    const Ptr<const T> end = begin + numberOfItems;
    for (; begin<end; ++begin)
        Serializer::deserializeImpl(r_stream, *begin);
    CIE_CHECK(!r_stream.fail(), "Stream in invalid state after array deserialization")
}


template <>
template <concepts::TriviallySerializable T>
inline void Serializer<tags::Binary>::serializeImpl(Ref<SerializerStream> r_stream, T object)
{
    r_stream.write(reinterpret_cast<const char*>(&object), sizeof(T));
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <concepts::NonTriviallySerializable T>
inline void Serializer<TTag>::serializeImpl(Ref<SerializerStream> r_stream, Ref<const T> r_object)
{
    r_object.serialize(r_stream, TTag());
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <concepts::Deserializable T>
inline void Serializer<TTag>::deserializeImpl(Ref<DeserializerStream> r_stream, Ref<T> r_output)
{
    if constexpr (concepts::TriviallyDeserializable<T>) {
        r_stream.read(reinterpret_cast<char*>(std::addressof(r_output)), sizeof(T));
    } else if constexpr (concepts::NonTriviallyDeserializable<T>) {
        T::deserialize(r_stream, r_output, TTag());
    }
}


} // namespace cie::io

#endif
