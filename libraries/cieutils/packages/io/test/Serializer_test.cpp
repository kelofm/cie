// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/Serializer.hpp"


namespace cie::io {


struct NonSerializableTestObject
{
    NonSerializableTestObject() {}

    NonSerializableTestObject(bool b, int i, long double d)
        : _bool(b), _int(i), _double(d)
    {}

    NonSerializableTestObject& operator=(int i)
    {
        _bool = bool(i);
        _int = i;
        _double = i*i;
        return *this;
    }

    friend bool operator==(Ref<const NonSerializableTestObject> r_left, Ref<const NonSerializableTestObject> r_right)
    {
        return (r_left._bool == r_right._bool) && (r_left._int == r_right._int) && (r_left._double == r_right._double);
    }

    friend std::ostream& operator<<(std::ostream& r_stream, Ref<const NonSerializableTestObject> r_object)
    {
        return r_stream << r_object._bool << ' ' << r_object._int << ' ' << r_object._double;
    }

    bool _bool;

    int _int;

    long double _double;
}; // class NonSerializableTestObject


struct TriviallySerializableTestObject : NonSerializableTestObject, io::TriviallySerializable
{
    using NonSerializableTestObject::NonSerializableTestObject;

    using NonSerializableTestObject::operator=;
};


struct SerializableTestObject : NonSerializableTestObject
{
    using NonSerializableTestObject::NonSerializableTestObject;

    using NonSerializableTestObject::operator=;

    void serialize(Ref<Traits::SerializerStream> r_stream,
                   tags::Binary = tags::Binary()) const
    {
        Serializer<tags::Binary> serializer;
        serializer.serialize(r_stream, _bool);
        serializer.serialize(r_stream, _int);
        serializer.serialize(r_stream, _double);
    }

    static void deserialize(Ref<Traits::DeserializerStream> r_stream,
                            Ref<NonSerializableTestObject> r_output,
                            tags::Binary = tags::Binary())
    {
        Serializer<tags::Binary> serializer;
        serializer.deserialize(r_stream, r_output._bool);
        serializer.deserialize(r_stream, r_output._int);
        serializer.deserialize(r_stream, r_output._double);
    }
};


template <class T>
void serializeDeserialize(Ref<const T> r_input, Ref<std::stringstream> r_stream)
{
    CIE_TEST_REQUIRE_NOTHROW(Serializer<tags::Binary>::serialize(r_stream, r_input));

    T deserialized;
    Serializer<tags::Binary>::deserialize<T>(r_stream, deserialized);
    CIE_TEST_CHECK(deserialized == r_input);
}


template <class T>
void serializeDeserialize(Ref<const T> r_input)
{
    std::stringstream stream;
    serializeDeserialize(r_input, stream);
}


CIE_TEST_CASE("Serializer", "[io]")
{
    CIE_TEST_CASE_INIT("Serializer")

    BinarySerializer serializer;

    #define CIE_TEST_SERIALIZER_FOR_TYPE(TYPE)                  \
        {                                                       \
            CIE_TEST_CASE_INIT(#TYPE)                           \
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
            serializer.deserialize(STREAM,BUFFER,64);           \
            for (int i=0; i<64; ++i)                            \
            {                                                   \
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
