// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/Serializer.hpp"


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

    static void deserialize(
        Ref<Traits::DeserializerStream> rStream,
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


} // namespace cie::io
