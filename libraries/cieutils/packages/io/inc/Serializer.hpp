#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/io/inc/ByteDeserializer.hpp"
#include "packages/io/inc/Traits.hpp"
#include "packages/types/inc/tags.hpp"
#include "packages/types/inc/types.hpp"


namespace cie::io {


/// @ingroup cieutils
template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
class Serializer {
public:
    using SerializerStream = Traits::SerializerStream;

    using DeserializerStream = Traits::DeserializerStream;

public:
    template <TriviallySerializable T>
    static void serialize(
        Ref<SerializerStream> rStream,
        T instance);

    template <NonTriviallySerializable T>
    static void serialize(
        Ref<SerializerStream> rStream,
        Ref<const T> rInstance);

    template <class T>
    requires Serializable<T,TTag>
    static void serialize(
        Ref<SerializerStream> rStream,
        Ptr<const T> begin,
        Size numberOfItems);

    template <class T, class TAllocator>
    requires Deserializable<T,TAllocator>
    static void deserialize(
        Ref<DeserializerStream> rStream,
        Ref<T> rOutput,
        Ref<TAllocator> rAllocator);

    template <class T, class TAllocator>
    requires Deserializable<T,TAllocator,TTag>
    static void deserialize(
        Ref<DeserializerStream> rStream,
        Ptr<T> begin,
        Size numberOfItems,
        Ref<TAllocator> rAllocator);

    /// @brief Deserialize from any stream satisfying @ref BinaryDeserializerStreamLike,
    ///        which includes @ref ByteDeserializer usable in SYCL kernels.
    template <class T, class TAllocator, class TStream>
    requires (
        std::is_same_v<TTag,tags::Binary>
        && BinaryDeserializerStreamLike<TStream>)
    static void deserialize(
        Ref<TStream> rStream,
        Ref<T> rOutput,
        Ref<TAllocator> rAllocator);

    /// @brief Deserialize an array of items from any stream satisfying
    ///        @ref BinaryDeserializerStreamLike, which includes @ref ByteDeserializer
    ///        usable in SYCL kernels.
    template <class T, class TAllocator, class TStream>
    requires (
        std::is_same_v<TTag,tags::Binary>
        && BinaryDeserializerStreamLike<TStream>)
    static void deserialize(
        Ref<TStream> rStream,
        Ptr<T> begin,
        Size numberOfItems,
        Ref<TAllocator> rAllocator);

private:
    template <TriviallySerializable T>
    static void serializeImpl(
        Ref<SerializerStream> rStream,
        T instance);

    template <NonTriviallySerializable T>
    static void serializeImpl(
        Ref<SerializerStream> rStream,
        Ref<const T> rInstance);

    template <class T, class TAllocator>
    requires Deserializable<T,TAllocator>
    static void deserializeImpl(
        Ref<DeserializerStream> rStream,
        Ref<T> rOutput,
        Ref<const TAllocator> rAllocator);
}; // class Serializer


/// @ingroup cieutils
using BinarySerializer = Serializer<tags::Binary>;


/// @ingroup cieutils
using TextSerializer = Serializer<tags::Text>;


} // namespace cie::io

#include "packages/io/impl/Serializer_impl.hpp"
