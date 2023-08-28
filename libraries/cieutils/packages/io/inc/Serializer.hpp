#ifndef CIE_UTILS_SERIALIZER_HPP
#define CIE_UTILS_SERIALIZER_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/io/inc/Traits.hpp"
#include "packages/types/inc/tags.hpp"
#include "packages/types/inc/types.hpp"

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
    static void serialize(Ref<SerializerStream> r_stream, T object);

    template <concepts::NonTriviallySerializable T>
    static void serialize(Ref<SerializerStream> r_stream, Ref<const T> r_object);

    template <class T>
    requires concepts::Serializable<T,TTag>
    static void serialize(Ref<SerializerStream> r_stream,
                          Ptr<const T> begin,
                          Size numberOfItems);

    template <concepts::Deserializable T>
    static void deserialize(Ref<DeserializerStream> r_stream, Ref<T> r_output);

    template <concepts::Deserializable<TTag> T>
    static void deserialize(Ref<DeserializerStream> r_stream,
                            Ptr<T> begin,
                            Size numberOfItems);

private:
    template <concepts::TriviallySerializable T>
    static void serializeImpl(Ref<SerializerStream> r_stream, T object);

    template <concepts::NonTriviallySerializable T>
    static void serializeImpl(Ref<SerializerStream> r_stream, Ref<const T> r_object);

    template <concepts::Deserializable T>
    static void deserializeImpl(Ref<DeserializerStream> r_stream, Ref<T> r_output);
}; // class Serializer


using BinarySerializer = Serializer<tags::Binary>;


using TextSerializer = Serializer<tags::Text>;


} // namespace cie::io


#include "packages/io/impl/Serializer_impl.hpp"

#endif
