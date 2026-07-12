#pragma once

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::io {


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <TriviallySerializable T>
void Serializer<TTag>::serialize(
    Ref<SerializerStream> rStream,
    T instance) {
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state before trivial serialization")
        #endif
        Serializer::serializeImpl(rStream, instance);
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state after trivial serialization")
        #endif
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <NonTriviallySerializable T>
void Serializer<TTag>::serialize(
    Ref<SerializerStream> rStream,
    Ref<const T> rInstance) {
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state before serialization")
        #endif
        Serializer::serializeImpl(rStream, rInstance);
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state after serialization")
        #endif
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <class T>
requires Serializable<T,TTag>
void Serializer<TTag>::serialize(
    Ref<SerializerStream> rStream,
    Ptr<const T> begin,
    Size numberOfItems) {
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state before array serialization")
        #endif
        const Ptr<const T> end = begin + numberOfItems;
        for (; begin<end; ++begin)
            Serializer::serialize(rStream, *begin);
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state after array serialization")
        #endif
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <class T, class TAllocator>
requires Deserializable<T,TAllocator>
void Serializer<TTag>::deserialize(
    Ref<DeserializerStream> rStream,
    Ref<T> rOutput,
    Ref<TAllocator> rAllocator) {
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state before trivial deserialization")
        #endif
        Serializer::deserializeImpl(
            rStream,
            rOutput,
            rAllocator);
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state after trivial deserialization")
        #endif
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <class T, class TAllocator>
requires Deserializable<T,TAllocator,TTag>
void Serializer<TTag>::deserialize(
    Ref<DeserializerStream> rStream,
    Ptr<T> begin,
    Size numberOfItems,
    Ref<TAllocator> rAllocator) {
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state before array deserialization")
        #endif
        const Ptr<const T> end = begin + numberOfItems;
        for (; begin<end; ++begin)
            Serializer::deserializeImpl(
                rStream,
                *begin,
                rAllocator);
        #ifndef NDEBUG
            CIE_CHECK(!rStream.fail(), "Stream in invalid state after array deserialization")
        #endif
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <TriviallySerializable T>
void Serializer<TTag>::serializeImpl(Ref<SerializerStream> rStream, T instance) {
    rStream.write(reinterpret_cast<const char*>(&instance), sizeof(T));
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <NonTriviallySerializable T>
void Serializer<TTag>::serializeImpl(Ref<SerializerStream> rStream, Ref<const T> rInstance) {
    rInstance.serialize(rStream, TTag());
}


template <concepts::AnyOf<tags::Binary,tags::Text> TTag>
template <class T, class TAllocator>
requires Deserializable<T,TAllocator>
void Serializer<TTag>::deserializeImpl(
    Ref<DeserializerStream> rStream,
    Ref<T> rOutput,
    Ref<TAllocator> rAllocator) {
    if constexpr (TriviallyDeserializable<T>) {
        rStream.read(reinterpret_cast<char*>(std::addressof(rOutput)), sizeof(T));
    } else if constexpr (NonTriviallyDeserializable<T,TAllocator>) {
        T::deserialize(
            rStream,
            rOutput,
            rAllocator,
            TTag());
    } else static_assert(TriviallyDeserializable<T>, "this type does not support deserialization with the provided allocator");
}


} // namespace cie::io
