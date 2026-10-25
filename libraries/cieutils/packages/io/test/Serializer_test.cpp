// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/ByteDeserializer.hpp"
#include "packages/io/inc/Serializer.hpp"

// --- STL Includes ---
#include <array>
#include <istream>


namespace cie::io {


struct NonSerializableTestObject {
    NonSerializableTestObject() {}

    NonSerializableTestObject(bool b, int i, long double d)
        : _bool(b), _int(i), _double(d)
    {}

    NonSerializableTestObject& operator=(int i) {
        _bool = bool(i);
        _int = i;
        _double = i*i;
        return *this;
    }

    friend bool operator==(
        Ref<const NonSerializableTestObject> rLeft,
        Ref<const NonSerializableTestObject> rRight) {
            return (rLeft._bool == rRight._bool) && (rLeft._int == rRight._int) && (rLeft._double == rRight._double);
    }

    friend std::ostream& operator<<(
        std::ostream& rStream,
        Ref<const NonSerializableTestObject> rObject) {
        return rStream << rObject._bool << ' ' << rObject._int << ' ' << rObject._double;
    }

    bool _bool;

    int _int;

    long double _double;
}; // class NonSerializableTestObject


struct TriviallySerializableTestObject : NonSerializableTestObject, io::TriviallySerializableBase {
    using NonSerializableTestObject::NonSerializableTestObject;

    using NonSerializableTestObject::operator=;
};


struct SerializableTestObject : NonSerializableTestObject {
    using NonSerializableTestObject::NonSerializableTestObject;

    using NonSerializableTestObject::operator=;

    void serialize(
        Ref<Traits::SerializerStream> rStream,
        tags::Binary = tags::Binary()) const {
            Serializer<tags::Binary> serializer;
            serializer.serialize(rStream, _bool);
            serializer.serialize(rStream, _int);
            serializer.serialize(rStream, _double);
    }

    /// @brief Deserialization from any stream type satisfying
    ///        @ref BinaryDeserializerStreamLike, including @ref ByteDeserializer
    ///        usable in SYCL kernels.
    template <class TStream>
    static void deserialize(
        Ref<TStream> rStream,
        Ref<NonSerializableTestObject> rOutput,
        std::allocator<std::byte> allocator,
        tags::Binary = tags::Binary()) {
            Serializer<tags::Binary> serializer;
            serializer.deserialize(
                rStream,
                rOutput._bool,
                allocator);
            serializer.deserialize(
                rStream,
                rOutput._int,
                allocator);
            serializer.deserialize(
                rStream,
                rOutput._double,
                allocator);
    }
};


template <class T>
void serializeDeserialize(Ref<const T> rInput, Ref<std::stringstream> rStream) {
    CIE_TEST_REQUIRE_NOTHROW(Serializer<tags::Binary>::serialize(rStream, rInput));

    std::allocator<int> allocator;

    T deserialized;
    Serializer<tags::Binary>::deserialize<T>(
        rStream,
        deserialized,
        allocator);
    CIE_TEST_CHECK(deserialized == rInput);
}


template <class T>
void serializeDeserialize(Ref<const T> rInput) {
    std::stringstream stream;
    serializeDeserialize(rInput, stream);
}


CIE_TEST_CASE("Serializer", "[io]")
{
    CIE_TEST_CASE_INIT("Serializer")

    BinarySerializer serializer;

    #define CIE_TEST_SERIALIZER_FOR_TYPE(TYPE)                  \
        {                                                       \
            CIE_TEST_CASE_INIT(#TYPE)                           \
            std::allocator<char> allocator;                     \
            TYPE VALUE;                                         \
            for (int i=-128; i<127; ++i)                        \
            {                                                   \
                VALUE = i;                                      \
                serializeDeserialize(VALUE);                    \
            }                                                   \
                                                                \
            TYPE BUFFER[64];                                    \
            for (int i=0; i<64; ++i)                            \
                BUFFER[i] = i;                                  \
                                                                \
            std::stringstream STREAM;                           \
            serializer.serialize(STREAM,BUFFER,64);             \
                                                                \
            for (int i=0; i<64; ++i)                            \
                BUFFER[i] = i + 128;                            \
                                                                \
            serializer.deserialize(                             \
                STREAM,                                         \
                BUFFER,                                         \
                64,                                             \
                allocator);                                     \
            for (int i=0; i<64; ++i) {                          \
                TYPE REFERENCE_VALUE;                           \
                REFERENCE_VALUE = i;                            \
                CIE_TEST_CHECK(BUFFER[i] == REFERENCE_VALUE);   \
            }                                                   \
        }

    CIE_TEST_SERIALIZER_FOR_TYPE(bool)
    CIE_TEST_SERIALIZER_FOR_TYPE(char)
    CIE_TEST_SERIALIZER_FOR_TYPE(int)
    CIE_TEST_SERIALIZER_FOR_TYPE(float)
    CIE_TEST_SERIALIZER_FOR_TYPE(double)
    CIE_TEST_SERIALIZER_FOR_TYPE(unsigned char)
    CIE_TEST_SERIALIZER_FOR_TYPE(unsigned int)
    CIE_TEST_SERIALIZER_FOR_TYPE(Size)
    CIE_TEST_SERIALIZER_FOR_TYPE(TriviallySerializableTestObject)
    CIE_TEST_SERIALIZER_FOR_TYPE(SerializableTestObject)
}


CIE_TEST_CASE("ByteDeserializer", "[io]")
{
    CIE_TEST_CASE_INIT("ByteDeserializer")

    static_assert(BinaryDeserializerStreamLike<std::istream>);
    static_assert(BinaryDeserializerStreamLike<std::stringstream>);
    static_assert(BinaryDeserializerStreamLike<ByteDeserializer>);

    std::allocator<std::byte> allocator;

    // Round-trip trivially serializable objects through a byte range.
    {
        std::stringstream stream;
        const TriviallySerializableTestObject input(false, 17, 2.71);
        Serializer<tags::Binary>::serialize(stream, input);
        const std::string bytes = stream.str();

        ByteDeserializer deserializer(
            reinterpret_cast<Ptr<const std::byte>>(bytes.data()),
            bytes.size());
        TriviallySerializableTestObject output;
        Serializer<tags::Binary>::deserialize(
            deserializer,
            output,
            allocator);
        CIE_TEST_CHECK(output == input);
        CIE_TEST_CHECK(deserializer.remaining() == 0ul);
        CIE_TEST_CHECK(deserializer.consumed() == bytes.size());
    }

    // Round-trip non-trivially serializable objects through a byte range.
    {
        std::stringstream stream;
        const SerializableTestObject input(true, -7, 3.14);
        Serializer<tags::Binary>::serialize(stream, input);
        const std::string bytes = stream.str();

        ByteDeserializer deserializer(
            reinterpret_cast<Ptr<const std::byte>>(bytes.data()),
            bytes.size());
        SerializableTestObject output;
        Serializer<tags::Binary>::deserialize(
            deserializer,
            output,
            allocator);
        CIE_TEST_CHECK(output == input);
        CIE_TEST_CHECK(deserializer.remaining() == 0ul);
    }

    // Deserialize an array of items from a byte range.
    {
        std::stringstream stream;
        std::array<int,8> input;
        for (std::size_t iItem=0ul; iItem<input.size(); ++iItem) input[iItem] = 3 * iItem;
        Serializer<tags::Binary>::serialize(
            stream,
            input.data(),
            input.size());
        const std::string bytes = stream.str();

        ByteDeserializer deserializer(
            reinterpret_cast<Ptr<const std::byte>>(bytes.data()),
            bytes.size());
        std::array<int,8> output;
        Serializer<tags::Binary>::deserialize(
            deserializer,
            output.data(),
            output.size(),
            allocator);
        for (std::size_t iItem=0ul; iItem<output.size(); ++iItem)
            CIE_TEST_CHECK(output[iItem] == input[iItem]);
    }
}


} // namespace cie::io
