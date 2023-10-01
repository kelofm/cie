#ifndef CIE_UTILS_SERIALIZER_HPP
#define CIE_UTILS_SERIALIZER_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/io/inc/Traits.hpp"
#include "packages/types/inc/tags.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <string>
#include <tuple>


namespace cie::io {


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
class Serializer
{
public:
    using SerializerStream = Traits::SerializerStream;

    using DeserializerStream = Traits::DeserializerStream;

public:
    template <concepts::TriviallySerializable T>
    static void serialize(Ref<SerializerStream> r_stream, T object)
    {
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state before trivial serialization")
        Serializer::serializeImpl(r_stream, object);
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state after trivial serialization")
    }

    template <concepts::NonTriviallySerializable T>
    static void serialize(Ref<SerializerStream> r_stream, Ref<const T> r_object)
    {
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state before serialization")
        Serializer::serializeImpl(r_stream, r_object);
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state after serialization")
    }

    template <class T>
    requires concepts::Serializable<T,TTag>
    static void serialize(Ref<SerializerStream> r_stream,
                          Ptr<const T> begin,
                          Size numberOfItems)
    {
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state before array serialization")
        const Ptr<const T> end = begin + numberOfItems;
        for (; begin<end; ++begin)
            Serializer::serialize(r_stream, *begin);
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state after array serialization")
    }

    template <concepts::Deserializable T>
    static void deserialize(Ref<DeserializerStream> r_stream, Ref<T> r_output)
    {
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state before trivial deserialization")
        Serializer::deserializeImpl(r_stream, r_output);
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state after trivial deserialization")
    }

    template <concepts::Deserializable<TTag> T>
    static void deserialize(Ref<DeserializerStream> r_stream,
                            Ptr<T> begin,
                            Size numberOfItems)
    {
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state before array deserialization")
        const Ptr<const T> end = begin + numberOfItems;
        for (; begin<end; ++begin)
            Serializer::deserializeImpl(r_stream, *begin);
        CIE_CHECK(!r_stream.fail(), "Stream in invalid state after array deserialization")
    }

private:
    template <concepts::TriviallySerializable T>
    static void serializeImpl(Ref<SerializerStream> r_stream, T object)
    {
        r_stream.write(reinterpret_cast<const char*>(&object), sizeof(T));
    }

    template <concepts::NonTriviallySerializable T>
    static void serializeImpl(Ref<SerializerStream> r_stream, Ref<const T> r_object)
    {
        r_object.serialize(r_stream, TTag());
    }

    template <concepts::Deserializable T>
    static void deserializeImpl(Ref<DeserializerStream> r_stream, Ref<T> r_output)
    {
        if constexpr (concepts::TriviallyDeserializable<T>) {
            r_stream.read(reinterpret_cast<char*>(std::addressof(r_output)), sizeof(T));
        } else if constexpr (concepts::NonTriviallyDeserializable<T>) {
            T::deserialize(r_stream, r_output, TTag());
        }
    }
}; // class Serializer


using BinarySerializer = Serializer<tags::Binary>;


using TextSerializer = Serializer<tags::Text>;


} // namespace cie::io

#endif
