// --- External Includes ---
#include "nlohmann/json.hpp"
#include "nlohmann/json-schema.hpp"

// --- Utiltity Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/io/inc/json.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <fstream>
#include <vector>
#include <array>
#include <iostream>
#include <format>


namespace cie::io {

// ----------------------------------------------------------
// ITERATOR
// ----------------------------------------------------------

namespace jsonimpl {
using iterator = nlohmann::json::iterator;
using const_iterator = nlohmann::json::const_iterator;

iterator advance(iterator it, int difference) {
    for (int i=0; i<difference; ++i)
        ++it;
    return it;
}

const_iterator advance(const_iterator it, int difference) {
    for (int i=0; i<difference; ++i)
        ++it;
    return it;
}
} // namespace detail


template <class Value>
JSONObject::IteratorBase<Value>::IteratorBase(value_type& rJSON)
    : _p_base(&rJSON),
        _index(0)
{}


template <class Value>
JSONObject::IteratorBase<Value>::IteratorBase(value_type& rJSON, difference_type index)
    : _p_base(&rJSON),
        _index(index)
{}


template <class T>
typename JSONObject::IteratorBase<T>::value_type*
JSONObject::IteratorBase<T>::get() {
    return _p_base;
}


template <class T>
const typename JSONObject::IteratorBase<T>::value_type*
JSONObject::IteratorBase<T>::get() const {
    return _p_base;
}


template <class Value>
typename JSONObject::IteratorBase<Value>::value_type
JSONObject::IteratorBase<Value>::operator*() {
    return JSONObject(&*jsonimpl::advance(this->get()->_pContents->begin(), _index));
}


template <class Value>
const typename JSONObject::IteratorBase<Value>::value_type
JSONObject::IteratorBase<Value>::operator*() const {
    return JSONObject(&*jsonimpl::advance(this->get()->_pContents->begin(), _index));
}


template <class Value>
JSONObject::IteratorBase<Value>&
JSONObject::IteratorBase<Value>::operator++() {
    ++_index;
    return *this;
}


template <class Value>
JSONObject::IteratorBase<Value>
JSONObject::IteratorBase<Value>::operator++(int)
{
    return JSONObject::IteratorBase<Value>(*this->get(), _index + 1);
}


template <class Value>
JSONObject::IteratorBase<Value>&
JSONObject::IteratorBase<Value>::operator--()
{
    --_index;
    return *this;
}


template <class Value>
JSONObject::IteratorBase<Value>
JSONObject::IteratorBase<Value>::operator--(int)
{
    return JSONObject::IteratorBase<Value>(*this->get(), _index - 1);
}


template <class Value>
bool
JSONObject::IteratorBase<Value>::operator!=(const IteratorBase& rRHS) const {
    return jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index)
    !=
    jsonimpl::advance(
        rRHS.get()->_pContents->begin(),
        rRHS._index);
}


template <class Value>
bool
JSONObject::IteratorBase<Value>::operator<(const IteratorBase& rRHS) const {
    return jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index)
    <
    jsonimpl::advance(
        rRHS.get()->_pContents->begin(),
        rRHS._index);
}


template <class Value>
bool
JSONObject::IteratorBase<Value>::operator<=(const IteratorBase& rRHS) const {
    return jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index)
    <=
    jsonimpl::advance(
        rRHS.get()->_pContents->begin(),
        rRHS._index);
}


template <class Value>
bool
JSONObject::IteratorBase<Value>::operator>(const IteratorBase& rRHS) const {
    return jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index)
    >
    jsonimpl::advance(
        rRHS.get()->_pContents->begin(),
        rRHS._index);
}


template <class Value>
bool
JSONObject::IteratorBase<Value>::operator>=(const IteratorBase& rRHS) const {
    return jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index)
    >=
    jsonimpl::advance(
        rRHS.get()->_pContents->begin(),
        rRHS._index);
}


template <class Value>
std::string
JSONObject::IteratorBase<Value>::key() const {
    auto it = jsonimpl::advance(
        this->get()->_pContents->begin(),
        _index);
    return it.key();
}


template <class Value>
typename JSONObject::IteratorBase<Value>::value_type
JSONObject::IteratorBase<Value>::value() {
    return **this;
}


template <class Value>
typename JSONObject::IteratorBase<Value>::value_type
JSONObject::IteratorBase<Value>::value() const {
    return **this;
}


template class JSONObject::IteratorBase<JSONObject>;


template class JSONObject::IteratorBase<const JSONObject>;

// ----------------------------------------------------------
// MEMBER TEMPLATE HELPERS
// ----------------------------------------------------------

template <class ValueType>
struct SetGetTemplate {
    static void set(
        nlohmann::json& rJSON,
        const ValueType& rValue) {
            CIE_BEGIN_EXCEPTION_TRACING
                rJSON = rValue;
            CIE_END_EXCEPTION_TRACING
    }

    static ValueType as(const nlohmann::json& rJSON) {
        CIE_BEGIN_EXCEPTION_TRACING
            return rJSON.get<ValueType>();
        CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct SetGetTemplate<JSONObject> {
    static void set(
        nlohmann::json& rJSON,
        const JSONObject& rValue) {
            CIE_BEGIN_EXCEPTION_TRACING
                rJSON = rValue.contents();
            CIE_END_EXCEPTION_TRACING
    }

    static JSONObject as(const nlohmann::json& rJSON) {
        CIE_BEGIN_EXCEPTION_TRACING
            return JSONObject(rJSON.dump());
        CIE_END_EXCEPTION_TRACING
    }
};


template <concepts::io::SupportedType ValueType>
struct TypeQueryTemplate {};


template <>
struct TypeQueryTemplate<bool> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_boolean();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = false;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<int> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_integer();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<long int> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_integer();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<long long> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_integer();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<Size> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_unsigned() && rJSON.contents().is_number_integer();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<float> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_float();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<double> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_number_float();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = 0;
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<std::string> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.contents().is_string();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
        CIE_BEGIN_EXCEPTION_TRACING
            rContents[rKey] = "";
            return JSONObject(&rContents[rKey]);
        CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<JSONObject> {
    static bool is(const JSONObject& rJSON) {
        return rJSON.isObject();
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = {};
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


template <>
struct TypeQueryTemplate<std::nullptr_t> {
    static bool is(const JSONObject& rJSON)
    { return rJSON.contents().is_null(); }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = JSONObject::content_type();
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


// ----------------------------------------------------------
// ARRAY MEMBER TEMPLATE HELPERS
// ----------------------------------------------------------

// Helper for array assignments
// WARNING: an adequately preallocated contiguous array is required
template <class ValueType>
struct AssignToContiguousArray {
    static void as(
        const nlohmann::json& rJSON,
        ValueType* pBegin) {
            CIE_BEGIN_EXCEPTION_TRACING
                for (const auto& rComponent : rJSON)
                    *pBegin++ = SetGetTemplate<ValueType>::as(rComponent);
            CIE_END_EXCEPTION_TRACING
    }
};


// Helper for std::vector<ValueType>
template <class ValueType>
struct SetGetTemplate<std::vector<ValueType>>
{
    static void set(nlohmann::json& rJSON,
                     const std::vector<ValueType>& rValue)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        rJSON.clear();

        for (const auto& rComponent : rValue)
            rJSON.push_back(rComponent);

        CIE_END_EXCEPTION_TRACING
    }

    static std::vector<ValueType> as(const nlohmann::json& rJSON)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        std::vector<ValueType> output;
        output.resize(rJSON.size());

        if (!output.empty())
            AssignToContiguousArray<ValueType>::as(rJSON, &output[0]);

        return output;

        CIE_END_EXCEPTION_TRACING
    }
};


// Helper for std::vector<JSONObject>
template <>
struct SetGetTemplate<std::vector<JSONObject>>
{
    static void set(nlohmann::json& rJSON,
                     const std::vector<JSONObject>& rValue)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        rJSON.clear();

        for (const auto& rComponent : rValue)
            rJSON.push_back(rComponent.contents());

        CIE_END_EXCEPTION_TRACING
    }

    static std::vector<JSONObject> as(const nlohmann::json& rJSON)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        std::vector<JSONObject> output;
        output.resize(rJSON.size());

        if (!output.empty())
            AssignToContiguousArray<JSONObject>::as(rJSON, &output[0]);

        return output;

        CIE_END_EXCEPTION_TRACING
    }
};


// Helper for StaticArray<ValueType,ArraySize>
template <class ValueType, Size ArraySize>
struct SetGetTemplate<StaticArray<ValueType,ArraySize>>
{
    static void set(nlohmann::json& rJSON,
                     const StaticArray<ValueType,ArraySize>& rValue)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        rJSON.clear();

        for (const auto& rComponent : rValue)
            rJSON.push_back(rComponent);

        CIE_END_EXCEPTION_TRACING
    }

    static StaticArray<ValueType,ArraySize> as(const nlohmann::json& rJSON)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        CIE_CHECK(
            rJSON.size() == ArraySize,
            "destination size (" + std::to_string(ArraySize) + ") does not match source size (" + std::to_string(rJSON.size()) + ")"
        )

        StaticArray<ValueType,ArraySize> output;
        AssignToContiguousArray<ValueType>::as(rJSON, &output[0]);

        return output;

        CIE_END_EXCEPTION_TRACING
    }
};


// Helper for StaticArray<JSONObject,ArraySize>
template <Size ArraySize>
struct SetGetTemplate<StaticArray<JSONObject,ArraySize>>
{
    static void set(
        nlohmann::json& rJSON,
        const StaticArray<JSONObject,ArraySize>& rValue) {
            CIE_BEGIN_EXCEPTION_TRACING
                rJSON.clear();
                for (const auto& rComponent : rValue)
                    rJSON.push_back(rComponent.contents());
            CIE_END_EXCEPTION_TRACING
    }

    static StaticArray<JSONObject,ArraySize> as(const nlohmann::json& rJSON) {
        CIE_BEGIN_EXCEPTION_TRACING
            CIE_CHECK(
                rJSON.size() == ArraySize,
                "destination size (" + std::to_string(ArraySize) + ") does not match source size (" + std::to_string(rJSON.size()) + ")")
            StaticArray<JSONObject,ArraySize> output;
            AssignToContiguousArray<JSONObject>::as(rJSON, &output[0]);
            return output;
        CIE_END_EXCEPTION_TRACING
    }
};


template <class ValueType>
struct TypeQueryTemplate<std::vector<ValueType>> {
    static bool is(const JSONObject& rJSON) {
        CIE_BEGIN_EXCEPTION_TRACING
            if (!rJSON.isArray())
                return false;
            const Size size = rJSON.size();
            for (Size i_component=0; i_component<size; ++i_component)
                if (!rJSON[i_component].is<ValueType>())
                    return false;
            return true;
        CIE_END_EXCEPTION_TRACING
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = JSONObject::content_type::array();
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


// Helper for StaticArray<ValueType,ArraySize>
template <class ValueType, Size ArraySize>
struct TypeQueryTemplate<StaticArray<ValueType,ArraySize>> {
    static bool is(const JSONObject& rJSON) {
        return TypeQueryTemplate<std::vector<ValueType>>::is(rJSON) && rJSON.size() == ArraySize;
    }

    static JSONObject addDefault(
        JSONObject&,
        JSONObject::content_type& rContents,
        Ref<const std::string> rKey) {
            CIE_BEGIN_EXCEPTION_TRACING
                rContents[rKey] = JSONObject::content_type::array();
                return JSONObject(&rContents[rKey]);
            CIE_END_EXCEPTION_TRACING
    }
};


// ----------------------------------------------------------
// MEMBER TEMPLATE INSTANTIATIONS
// ----------------------------------------------------------

template <class ValueType>
void
JSONObject::SetGet<ValueType>::set(
    JSONObject& rJSON,
    const ValueType& rValue) {
        SetGetTemplate<ValueType>::set(
            *rJSON._pContents,
            rValue);
}


template <class ValueType>
ValueType
JSONObject::SetGet<ValueType>::as(const JSONObject& rJSON)
{
    return SetGetTemplate<ValueType>::as(*rJSON._pContents);
}


template <class ValueType>
bool JSONObject::TypeQuery<ValueType>::is(const JSONObject& rJSON) {
    return TypeQueryTemplate<ValueType>::is(rJSON);
}


template <class ValueType>
JSONObject JSONObject::TypeQuery<ValueType>::addDefault(JSONObject& rJSON,
                                                         Ref<const std::string> rKey)
{
    return TypeQueryTemplate<ValueType>::addDefault(
        rJSON,
        *rJSON._pContents.get(),
        rKey
    );
}


// Template instantiations
#define CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(ValueType)        \
    template struct JSONObject::SetGet<ValueType>;                      \
    template struct JSONObject::SetGet<std::vector<ValueType>>;         \
    template struct JSONObject::SetGet<StaticArray<ValueType,1>>;       \
    template struct JSONObject::SetGet<StaticArray<ValueType,2>>;       \
    template struct JSONObject::SetGet<StaticArray<ValueType,3>>;       \
    template struct JSONObject::TypeQuery<ValueType>;                   \
    template struct JSONObject::TypeQuery<std::vector<ValueType>>;      \
    template struct JSONObject::TypeQuery<StaticArray<ValueType,1>>;    \
    template struct JSONObject::TypeQuery<StaticArray<ValueType,2>>;    \
    template struct JSONObject::TypeQuery<StaticArray<ValueType,3>>;


CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(Size)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(int)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(long int)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(long long)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(float)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(double)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(std::string)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(JSONObject)
CIE_INSTANTIATE_JSON_TEMPLATE_SPECIALIZATIONS(std::nullptr_t)


// bool requires special treatment because of std::vector<bool> ¯\_(ツ)_/¯
template struct JSONObject::SetGet<bool>;
template struct JSONObject::SetGet<StaticArray<bool,1>>;
template struct JSONObject::SetGet<StaticArray<bool,2>>;
template struct JSONObject::SetGet<StaticArray<bool,3>>;
template struct JSONObject::TypeQuery<bool>;
template struct JSONObject::TypeQuery<StaticArray<bool,1>>;
template struct JSONObject::TypeQuery<StaticArray<bool,2>>;
template struct JSONObject::TypeQuery<StaticArray<bool,3>>;


// ----------------------------------------------------------
// CONSTRUCTOR / DESTRUCTOR
// ----------------------------------------------------------

JSONObject::JSONObject()
    :   _pContents(new nlohmann::json),
        _isRoot(true)
{}


JSONObject::JSONObject(const JSONObject& rRHS)
    :   _pContents(new nlohmann::json(*rRHS._pContents)),
        _isRoot(true)
{}


JSONObject::JSONObject(JSONObject&& rRHS)
    :   _pContents(rRHS._pContents),
        _isRoot(rRHS._isRoot) {
            rRHS._pContents = nullptr;
            rRHS._isRoot = false;
}


JSONObject& JSONObject::operator=(JSONObject&& rRHS) {
    if (&rRHS != this) {
        this->~JSONObject();
        _pContents = rRHS._pContents;
        _isRoot = rRHS._isRoot;
        rRHS._pContents = nullptr;
        rRHS._isRoot = false;
    }
    return *this;
}


JSONObject& JSONObject::operator=(const JSONObject& rRHS) {
    if (&rRHS == this) return *this;

    CIE_BEGIN_EXCEPTION_TRACING
        this->~JSONObject();
        _pContents = new JSONObject::content_type(*rRHS._pContents);
        _isRoot = true;
        return *this;
    CIE_END_EXCEPTION_TRACING
}


JSONObject::JSONObject(Ref<const std::string> rJSONString)
    :   _pContents(new JSONObject::content_type(JSONObject::content_type::parse(rJSONString))),
        _isRoot(true)
{}


JSONObject::JSONObject(std::string&& rJSONString)
    :   _pContents(new JSONObject::content_type(JSONObject::content_type::parse(std::move(rJSONString)))),
        _isRoot(true)
{}


JSONObject::JSONObject(const std::filesystem::path& rFilePath)
    : JSONObject() {
        CIE_BEGIN_EXCEPTION_TRACING
            std::ifstream file(rFilePath);
            try {
                file >> *_pContents;
            } catch (std::exception& rException) {
                CIE_THROW(
                    Exception,
                    std::format(
                        "Failed to parse '{}'\n{}",
                        rFilePath.string(),
                        rException.what()))
            }
        CIE_END_EXCEPTION_TRACING
}


JSONObject::JSONObject(std::istream& rStream)
    : JSONObject() {
        CIE_BEGIN_EXCEPTION_TRACING
            rStream >> *_pContents;
        CIE_END_EXCEPTION_TRACING
}


JSONObject::JSONObject(JSONObject::content_type* pContents)
        :   _pContents(pContents),
            _isRoot(false)
{}


JSONObject::JSONObject(const JSONObject::content_type* pContents)
        :   _pContents(pContents),
            _isRoot(false)
{}


JSONObject::~JSONObject() {
    if (_isRoot && _pContents) {
        delete _pContents;
        _pContents = nullptr;
        _isRoot = false;
    }
}


// ----------------------------------------------------------
// NON-TEMPLATE MEMBERS
// ----------------------------------------------------------

JSONObject JSONObject::operator[](Ref<const std::string> rKey) {
    CIE_CHECK(
        this->hasKey(rKey),
        "'" + rKey + "' is not a key")
    CIE_BEGIN_EXCEPTION_TRACING
        return JSONObject(&_pContents->operator[](rKey));
    CIE_END_EXCEPTION_TRACING
}


const JSONObject JSONObject::operator[](Ref<const std::string> rKey) const {
    CIE_CHECK(
        this->hasKey(rKey),
        "'" + rKey + "' is not a key")
    CIE_BEGIN_EXCEPTION_TRACING
        return JSONObject(&_pContents->operator[](rKey));
    CIE_END_EXCEPTION_TRACING
}


JSONObject JSONObject::operator[](Size index) {
    CIE_BEGIN_EXCEPTION_TRACING
        return JSONObject(&_pContents->operator[](index));
    CIE_END_EXCEPTION_TRACING
}


const JSONObject JSONObject::operator[](Size index) const {
    CIE_BEGIN_EXCEPTION_TRACING
        return JSONObject(&_pContents->operator[](index));
    CIE_END_EXCEPTION_TRACING
}


bool JSONObject::hasKey(Ref<const std::string> rKey) const {
    return _pContents->contains(rKey);
}


Size JSONObject::size() const {
    return _pContents->size();
}


JSONObject::iterator JSONObject::begin() {
    return JSONObject::iterator(*this);
}


JSONObject::const_iterator JSONObject::begin() const {
    return JSONObject::const_iterator(*this);
}


JSONObject::iterator JSONObject::end() {
    return JSONObject::iterator(*this, this->size());
}


JSONObject::const_iterator JSONObject::end() const {
    return JSONObject::const_iterator(*this, this->size());
}


bool JSONObject::isArray() const {
    return _pContents->is_array();
}


bool JSONObject::isObject() const {
    return _pContents->is_object();
}


const JSONObject::content_type& JSONObject::contents() const {
    return *_pContents;
}


JSONObject::content_type& JSONObject::contents() {
    return *_pContents;
}


JSONObject& JSONObject::add(
    Ref<const std::string> rKey,
    Ref<const std::string> rValue,
    bool allowOverwrite) {
        CIE_BEGIN_EXCEPTION_TRACING
            if (this->hasKey(rKey)) {
                if (!allowOverwrite)
                    CIE_THROW(Exception, "'" + rKey + "' is an existing key")
                this->operator[](rKey).set(rValue);
            } else {
                this->addDefault<std::string>(rKey).set(rValue);
            }
            return *this;
        CIE_END_EXCEPTION_TRACING
}


JSONObject& JSONObject::set(Ref<const std::string> rValue) {
    CIE_BEGIN_EXCEPTION_TRACING

    JSONObject::SetGet<std::string>::set(*this, rValue);
    return *this;

    CIE_END_EXCEPTION_TRACING
}


void JSONObject::prettyPrint(
    Ref<std::ostream> rStream,
    int indentation) const {
        rStream << this->contents().dump(indentation);
}


// ----------------------------------------------------------
// JSONSchema
// ----------------------------------------------------------


nlohmann::json_schema::schema_loader makeJSONSchemaLoader(Ref<const JSONSchemaLoader> rLoader) {
    return [&rLoader] (
        const nlohmann::json_uri& rURI,
        nlohmann::json& rResource) {
            CIE_BEGIN_EXCEPTION_TRACING
                JSONObject loadedResource = rLoader.load(rURI.path());
                rResource = std::move(loadedResource.contents());
            CIE_END_EXCEPTION_TRACING
        };
}


struct JSONSchema::Impl {
    JSONObject::content_type schema;

    nlohmann::json_schema::json_validator validator;

    JSONSchemaLoader loader;
}; // struct JSONSchema::Impl


JSONSchema::~JSONSchema() = default;


JSONSchema::JSONSchema()
    : _pImpl(new Impl)
{}


JSONSchema::JSONSchema(
    Ref<const JSONObject> rSchema,
    Ref<const JSONSchemaLoader> rLoader)
    : JSONSchema(JSONObject(rSchema), JSONSchemaLoader(rLoader))
{}


JSONSchema::JSONSchema(
    RightRef<JSONObject> rSchema,
    RightRef<JSONSchemaLoader> rLoader)
        : _pImpl(new Impl) {
            CIE_BEGIN_EXCEPTION_TRACING
                _pImpl->loader = std::move(rLoader);
                _pImpl->validator = nlohmann::json_schema::json_validator(makeJSONSchemaLoader(_pImpl->loader));
                _pImpl->schema = std::move(rSchema.contents());
                _pImpl->validator.set_root_schema(_pImpl->schema);
            CIE_END_EXCEPTION_TRACING
}


JSONSchema& JSONSchema::operator=(JSONSchema&&) noexcept = default;


enum class JSONSchemaConstraint {
    Property,
    AnyOf,
    OneOf,
    AllOf
}; // enum class JSONSchemaConstraint


class RecursiveFillErrorHandler : public nlohmann::json_schema::error_handler {
public:
    RecursiveFillErrorHandler()
        : _triggered(false)
    {}

    void error(
        const nlohmann::json::json_pointer&,
        const nlohmann::json&,
        Ref<const std::string>) override {
            _triggered = true;
    }

    [[nodiscard]] bool isTriggered() const noexcept {
        return _triggered;
    }

private:
    bool _triggered;
}; // RecursiveFillErrorHandler


bool recursiveFill(
    Ref<const JSONObject::content_type> rSchema,
    Ref<JSONObject::content_type> rJSON,
    Ref<const JSONSchemaLoader> rLoader,
    Ref<std::unordered_map<std::string,JSONObject>> rCache);


template <JSONSchemaConstraint TConstraint>
bool recursiveFillConstraint(
    Ref<const JSONObject::content_type> rSchema,
    Ref<JSONObject::content_type> rJSON,
    Ref<const JSONSchemaLoader> rLoader,
    Ref<std::unordered_map<std::string,JSONObject>> rCache) {
        bool modified = false;

        if constexpr (TConstraint == JSONSchemaConstraint::Property) {
            for (auto it=rSchema.begin(); it!=rSchema.end(); ++it) {
                // Handle references.
                std::optional<JSONObject> maybeLoadedSchema;
                if (it->contains("$ref")) {
                    Ref<const std::string> rURI = it->find("$ref")->get<std::string>();
                    const auto itCache = rCache.find(rURI);
                    if (itCache == rCache.end())
                        CIE_BEGIN_EXCEPTION_TRACING
                            maybeLoadedSchema = rCache.emplace(
                                rURI,
                                rLoader.load(rURI)
                            ).first->second;
                        CIE_END_EXCEPTION_TRACING
                    else
                        maybeLoadedSchema = itCache->second;
                }
                const auto& rResolvedSchema = maybeLoadedSchema.has_value()
                    ? maybeLoadedSchema.value().contents()
                    : *it;

                if (rResolvedSchema.contains("default")) {
                    const auto& rDefaultProperty = *rResolvedSchema.find("default");
                    if (!rJSON.contains(it.key())) {
                        modified = true;
                        rJSON.operator[](it.key()) = rDefaultProperty;
                    } else {
                        auto& rProperty = rJSON[it.key()];
                        if (rProperty.is_object()) {
                            modified |= recursiveFill(*it, rProperty, rLoader, rCache);
                        } else if (rProperty.is_array()) {
                            for (auto& rItem : rProperty)
                                modified |= recursiveFill(*it, rItem, rLoader, rCache);
                        }
                    }
                } /*property has default*/ else if (rJSON.contains(it.key())) {
                    auto& rProperty = rJSON[it.key()];
                    if (rProperty.is_structured()) {
                        modified |= recursiveFill(*it, rProperty, rLoader, rCache);
                    }
                } // json contains property
            } // for it
        } /*TConstraint == JSONSchemaConstraint::Property*/ else if constexpr (TConstraint == JSONSchemaConstraint::AnyOf) {
            // Loop over all subschemas, find out which ones apply to the
            // subject, and apply defaults to it if necessary.
            for (const auto& rSubschema : rSchema) {
                // Handle references.
                std::optional<JSONObject> maybeLoadedSchema;
                if (rSubschema.contains("$ref")) {
                    Ref<const std::string> rURI = rSubschema.find("$ref")->get<std::string>();
                    const auto itCache = rCache.find(rURI);
                    if (itCache == rCache.end())
                        CIE_BEGIN_EXCEPTION_TRACING
                            maybeLoadedSchema = rCache.emplace(
                                rURI,
                                rLoader.load(rURI)
                            ).first->second;
                        CIE_END_EXCEPTION_TRACING
                    else
                        maybeLoadedSchema = itCache->second;
                }
                const auto& rResolvedSchema = maybeLoadedSchema.has_value()
                    ? maybeLoadedSchema.value().contents()
                    : rSubschema;

                RecursiveFillErrorHandler errorHandler;
                nlohmann::json_schema::json_validator validator(makeJSONSchemaLoader(rLoader));
                validator.set_root_schema(rResolvedSchema);
                validator.validate(
                    rJSON,
                    errorHandler);
                if (!errorHandler.isTriggered()) modified |= recursiveFill(rResolvedSchema, rJSON, rLoader, rCache);
            } // for rSubschema
        } /*TConstraint == JSONSchemaConstraint::AnyOf*/ else if constexpr (TConstraint == JSONSchemaConstraint::OneOf) {
            // Loop over all subschemas, find out which one applies to
            // the subject, and apply defaults to it if necessary.
            for (const auto& rSubschema : rSchema) {
                // Handle references.
                std::optional<JSONObject> maybeLoadedSchema;
                if (rSubschema.contains("$ref")) {
                    Ref<const std::string> rURI = rSubschema.find("$ref")->get<std::string>();
                    const auto itCache = rCache.find(rURI);
                    if (itCache == rCache.end())
                        CIE_BEGIN_EXCEPTION_TRACING
                            maybeLoadedSchema = rCache.emplace(
                                rURI,
                                rLoader.load(rURI)
                            ).first->second;
                        CIE_END_EXCEPTION_TRACING
                    else
                        maybeLoadedSchema = itCache->second;
                }
                const auto& rResolvedSchema = maybeLoadedSchema.has_value()
                    ? maybeLoadedSchema.value().contents()
                    : rSubschema;

                RecursiveFillErrorHandler errorHandler;
                nlohmann::json_schema::json_validator validator(makeJSONSchemaLoader(rLoader));
                validator.set_root_schema(rResolvedSchema);
                validator.validate(
                    rJSON,
                    errorHandler);
                if (!errorHandler.isTriggered()) {
                    modified |= recursiveFill(rResolvedSchema, rJSON, rLoader, rCache);
                    break;
                }
            } // for rSubschema
        } /*TConstraint == JSONSchemaConstraint::OneOf*/ else if constexpr (TConstraint == JSONSchemaConstraint::AllOf) {
            for (const auto& rSubschema : rSchema) {
                // Handle references.
                std::optional<JSONObject> maybeLoadedSchema;
                if (rSubschema.contains("$ref")) {
                    Ref<const std::string> rURI = rSubschema.find("$ref")->get<std::string>();
                    const auto itCache = rCache.find(rURI);
                    if (itCache == rCache.end())
                        CIE_BEGIN_EXCEPTION_TRACING
                            maybeLoadedSchema = rCache.emplace(
                                rURI,
                                rLoader.load(rURI)
                            ).first->second;
                        CIE_END_EXCEPTION_TRACING
                    else
                        maybeLoadedSchema = itCache->second;
                }
                const auto& rResolvedSchema = maybeLoadedSchema.has_value()
                    ? maybeLoadedSchema.value().contents()
                    : rSubschema;

                modified |= recursiveFill(rResolvedSchema, rJSON, rLoader, rCache);
            } // for rSubschema
        } /*TConstraint == JSONSchemaConstraint::AllOf*/ else static_assert(TConstraint == JSONSchemaConstraint::Property, "invalid JSON constraint");

        return modified;
}


bool recursiveFill(
    Ref<const JSONObject::content_type> rSchema,
    Ref<JSONObject::content_type> rJSON,
    Ref<const JSONSchemaLoader> rLoader,
    Ref<std::unordered_map<std::string,JSONObject>> rCache) {
        bool modified = false;
        for (auto it=rSchema.begin(); it!=rSchema.end(); ++it) {
            Ref<const std::string> rConstraintName = it.key();
            if (
                   rConstraintName == "title"
                || rConstraintName == "description"
                || rConstraintName == "$id"
                || rConstraintName == "$schema"
                || rConstraintName == "additionalProperties"
                || rConstraintName == "required") {
                    // Not an actual constraint. Do nothing.
            } else if (rConstraintName == "type") {
                // Not a constraint default assignment has anything to do with.
            } else if (rConstraintName == "default") {
                // Handled elsewhere.
            } else if (rConstraintName == "properties") {
                modified |= recursiveFillConstraint<JSONSchemaConstraint::Property>(
                    *it,
                    rJSON,
                    rLoader,
                    rCache);
            } else if (rConstraintName == "oneOf") {
                modified |= recursiveFillConstraint<JSONSchemaConstraint::OneOf>(
                    *it,
                    rJSON,
                    rLoader,
                    rCache);
            } else if (rConstraintName == "anyOf") {
                modified |= recursiveFillConstraint<JSONSchemaConstraint::AnyOf>(
                    *it,
                    rJSON,
                    rLoader,
                    rCache);
            } else if (rConstraintName == "allOf") {
                modified |= recursiveFillConstraint<JSONSchemaConstraint::AllOf>(
                    *it,
                    rJSON,
                    rLoader,
                    rCache);
            } else if (rConstraintName == "items") {
                for (auto& rItem : rJSON)
                    modified |= recursiveFill(*it, rItem, rLoader, rCache);
            } else if (rConstraintName == "$ref") {
                Ref<const std::string> rURI = it->get<std::string>();
                auto itCache = rCache.find(rURI);
                if (itCache == rCache.end())
                    CIE_BEGIN_EXCEPTION_TRACING
                        itCache = rCache.emplace(
                            rURI,
                            rLoader.load(rURI)).first;
                    CIE_END_EXCEPTION_TRACING
                recursiveFill(itCache->second.contents(), rJSON, rLoader, rCache);
            } else CIE_THROW(Exception, std::format(
                "JSON constraint '{}' is not supported",
                rConstraintName))
        } // for it
        return modified;
}


void JSONSchema::validate(Ref<const JSONObject> rJSON) const {
    CIE_BEGIN_EXCEPTION_TRACING
        _pImpl->validator.validate(rJSON.contents());
    CIE_END_EXCEPTION_TRACING
}


void JSONSchema::validateAndFillDefaults(Ref<JSONObject> rJSON) const {
    CIE_BEGIN_EXCEPTION_TRACING
        std::unordered_map<std::string,JSONObject> cache;
        do {
            _pImpl->validator.validate(rJSON.contents());
        } while (recursiveFill(_pImpl->schema, rJSON.contents(), _pImpl->loader, cache));
    CIE_END_EXCEPTION_TRACING
}


void JSONSchema::resolve(Ref<const JSONSchemaLoader>) {
    CIE_THROW(NotImplementedException, "")
}


// ----------------------------------------------------------
// JSONSchema
// ----------------------------------------------------------


JSONSchemaLoader::JSONSchemaLoader()
    : JSONSchemaLoader(cie::getDataPath() / "data" / "schemas")
{}


JSONSchemaLoader::JSONSchemaLoader(Ref<const std::filesystem::path> rLibraryRoot)
    : JSONSchemaLoader(std::span<const std::filesystem::path>(&rLibraryRoot, 1))
{}


JSONSchemaLoader::JSONSchemaLoader(std::span<const std::filesystem::path> libraryRoots)
    : _roots() {
        CIE_BEGIN_EXCEPTION_TRACING
            for (Ref<const std::filesystem::path> rRoot : libraryRoots)
                this->addLibraryRoot(rRoot);
        CIE_END_EXCEPTION_TRACING
}


JSONSchemaLoader::~JSONSchemaLoader() = default;


void JSONSchemaLoader::addLibraryRoot(Ref<const std::filesystem::path> rLibraryRoot) {
    CIE_BEGIN_EXCEPTION_TRACING
        //this->validateLibraryRoot(rLibraryRoot);
        if (std::find(_roots.begin(), _roots.end(), rLibraryRoot) == _roots.end())
            _roots.push_back(rLibraryRoot);
    CIE_END_EXCEPTION_TRACING
}


JSONObject JSONSchemaLoader::load(Ref<const std::string> rURI) const {
    const std::filesystem::path relativeURIPath = rURI.empty()
        ? ""
        : rURI[0] == '/'
            ? rURI.substr(1, rURI.npos) + ".json"
            : rURI + ".json";

    CIE_BEGIN_EXCEPTION_TRACING
        std::optional<std::filesystem::path> maybePath;
        for (Ref<const std::filesystem::path> rRoot : _roots) {
            const std::filesystem::path candidate = rRoot / relativeURIPath;
            const auto status = std::filesystem::status(candidate);
            if (status.type() == std::filesystem::file_type::regular) {
                maybePath = candidate;
                break;
            } // @todo: deal with symlinks
        } // for rRoot in roots

        // Make sure it was found and loaded.
        CIE_CHECK(
            maybePath.has_value(),
            std::format(
                "'{}' does not refer to a schema in the following schema libraries:\n{}",
                rURI,
                [this](){
                    std::stringstream message;
                    for (Ref<const std::filesystem::path> rRoot : _roots)
                        message << "  " << rRoot << "\n";
                    return message.str();
                }()))

        return JSONObject(maybePath.value());
    CIE_END_EXCEPTION_TRACING
}


//void JSONSchemaLoader::validateLibraryRoot(Ref<const std::filesystem::path> rRoot) {
//    const auto status = std::filesystem::status(rRoot);
//        switch (status.type()) {
//            case std::filesystem::file_type::directory: break; // <== ok
//            case std::filesystem::file_type::regular:
//                CIE_THROW(
//                    Exception,
//                    std::format(
//                        "expecting a directory, but '{}' is a file",
//                        rRoot.string()))
//            case std::filesystem::file_type::none:
//                CIE_THROW(
//                    Exception,
//                    std::format(
//                        "expecting a directory, but '{}' does not exist",
//                        rRoot.string()))
//            case std::filesystem::file_type::symlink:
//                CIE_THROW(
//                    NotImplementedException,
//                    std::format(
//                        "expecting a directory, but '{}' is a symlink, which is not supported yet",
//                        rRoot.string()))
//            default:
//                CIE_THROW(
//                    Exception,
//                    std::format(
//                        "expecting a directory, but '{}' is of unsupported type",
//                        rRoot.string()))
//        } // switch status.type()
//}


// ----------------------------------------------------------
// NON-MEMBERS
// ----------------------------------------------------------

std::ostream& operator<<(
    std::ostream& rStream,
    const JSONObject& rJSON) {
        CIE_BEGIN_EXCEPTION_TRACING
            return rStream << rJSON.contents();
        CIE_END_EXCEPTION_TRACING
}


std::ostream& operator<<(
    [[maybe_unused]] std::ostream& rStream,
    [[maybe_unused]] Ref<const JSONSchema> rSchema) {
        CIE_THROW(NotImplementedException, "")
}


} // namespace cie::io