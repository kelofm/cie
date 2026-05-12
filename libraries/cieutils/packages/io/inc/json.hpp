#pragma once

// --- External Includes ---
#include "nlohmann/json_fwd.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/RuntimeConst.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/types/inc/modifiers.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <istream>
#include <ostream>
#include <filesystem>
#include <span>


namespace cie::io {
class JSONObject;
} // namespace cie::io


namespace cie::concepts::io {

template <class T>
concept SupportedBase
=  std::same_as<std::decay_t<T>, bool>
|| std::same_as<std::decay_t<T>, int>
|| std::same_as<std::decay_t<T>, long int>
|| std::same_as<std::decay_t<T>, long long>
|| std::same_as<std::decay_t<T>, Size>
|| std::same_as<std::decay_t<T>, float>
|| std::same_as<std::decay_t<T>, double>
|| std::same_as<std::decay_t<T>, std::string>
|| std::same_as<std::decay_t<T>, cie::io::JSONObject>;

template <class T>
concept SupportedType
=  SupportedBase<T>
|| (Container<std::decay_t<T>> && SupportedBase<typename std::decay_t<T>::value_type>);
} // namespace cie::concepts::io


namespace cie::io {


/** @brief Basic interface for a json parser
 *  @details supported types for io:
 *           - bool
 *           - int
 *           - cie::Size
 *           - float
 *           - double
 *           - std::string
 *           - JSONObject
 *           - std::vector of any above except bool
 *           - StaticArray of any above with size 1, 2, or 3
 *  @ingroup cieutils
 */
class JSONObject {
public:
    using content_type = nlohmann::json;

private:
    template <class Value>
    class IteratorBase {
    public:
        using value_type      = typename CopyConstQualifier<Value,JSONObject>::Type;
        using pointer         = value_type*;
        using reference       = value_type&;
        using difference_type = int;

    public:
        IteratorBase(value_type& rJSON);

        IteratorBase(value_type& rJSON, difference_type index);

        IteratorBase() = delete;

        IteratorBase(IteratorBase&& rRhs) = default;

        IteratorBase(const IteratorBase& rRhs)
        requires concepts::Const<Value> = default;

        IteratorBase& operator=(IteratorBase&& rRhs) = default;

        IteratorBase& operator=(const IteratorBase& rRhs)
        requires concepts::Const<Value> = default;

        value_type operator*();

        const value_type operator*() const;

        IteratorBase& operator++();

        IteratorBase operator++(int);

        IteratorBase& operator--();

        IteratorBase operator--(int);

        //IteratorBase& operator+=(difference_type difference);

        //IteratorBase& operator-=(difference_type difference);

        //IteratorBase operator+(difference_type rhs);

        //IteratorBase operator-(difference_type rhs);

        //bool operator==(const IteratorBase& rRhs);

        bool operator!=(const IteratorBase& rRhs);

        bool operator<(const IteratorBase& rRhs);

        bool operator<=(const IteratorBase& rRhs);

        bool operator>(const IteratorBase& rRhs);

        bool operator>=(const IteratorBase& rRhs);

        std::string key() const;

        value_type value();

        const value_type value() const;

    private:
        value_type*     _p_base;
        difference_type _index;
    }; // class IteratorBase

public:
    using iterator       = IteratorBase<JSONObject>;
    using const_iterator = IteratorBase<const JSONObject>;

public:
    JSONObject();

    JSONObject(const JSONObject& rRhs);

    JSONObject(JSONObject&& rRhs);

    JSONObject& operator=(JSONObject&& rRhs);

    JSONObject& operator=(const JSONObject& rRhs);

    JSONObject(const std::string& rJSONString);

    explicit JSONObject(std::string&& rJSONString);

    explicit JSONObject(const std::filesystem::path& r_filePath);

    JSONObject(std::istream& r_stream);

    /// Constructor for operator[]
    JSONObject(
        content_type* p_contents,
        const JSONObject* p_root);

    /// Constructor for operator[] const
    JSONObject(
        const content_type* p_contents,
        const JSONObject* p_root);

    ~JSONObject();

    /// Get the value associated to the specified key
    JSONObject operator[](const std::string& rKey);

    /// Get the value associated to the specified key
    const JSONObject operator[](const std::string& rKey) const;

    /// Get the specified component if this is an array
    JSONObject operator[](Size index);

    /// Get the specified component if this is an array
    const JSONObject operator[](Size index) const;

    /** Create a new string item with the specified key and value
     *  @details implicitly converts char arrays to std::string
     */
    JSONObject& add(
        const std::string& rKey,
        const std::string& r_value,
        bool allowOverwrite = false);

    /// Create a new item with the specified key and value
    template <concepts::io::SupportedType ValueType>
    JSONObject& add(
        const std::string& rKey,
        const ValueType& r_value,
        bool allowOverwrite = false);

    /** Set the json value to the specified string
     *  @details implicitly converts char arrays to std::string
     */
    JSONObject& set(const std::string& r_value);

    /// Set the json value to the specified one
    template <concepts::io::SupportedType ValueType>
    JSONObject& set(const ValueType& r_value);

    /// Parse this as an instance of the template parameter
    template <concepts::io::SupportedType ValueType>
    ValueType as() const;

    /// Return true if this can be converted to an instance of the template parameter
    template <concepts::io::SupportedType ValueType>
    bool is() const;

    /// Check if this is a json array
    bool isArray() const;

    /// Check if this is a json object
    bool isObject() const;

    /// Return true if this has an item with a matching key
    bool hasKey(const std::string& rKey) const;

    /// Item size
    Size size() const;

    iterator begin();

    const_iterator begin() const;

    iterator end();

    const_iterator end() const;

    /// Get wrapped object
    const content_type& contents() const;

    /// Get wrapped object
    content_type& contents();

    /// Get the root that holds the resources
    const JSONObject& root() const;

    void prettyPrint(
        Ref<std::ostream> rStream,
        int indentation = 4) const;

private:
    /// Set / get helper class
    template <class ValueType>
    struct SetGet {
        static void set(
            JSONObject& rJSON,
            const ValueType& r_value);

        static ValueType as(const JSONObject& rJSON);
    };

    /// Type identification and initialization helper
    template <class ValueType>
    struct TypeQuery {
        static bool is(const JSONObject& rJSON);

        static JSONObject addDefault(
            JSONObject& rJSON,
            const std::string& rKey);
    };

    /// Type-safe initialization of a new item's value
    template <concepts::io::SupportedType ValueType>
    JSONObject addDefault(const std::string& rKey);

protected:
    utils::RuntimeConst<content_type> _pContents;

private:
    const JSONObject* _pRoot;
}; // class JSONObject


/// @ingroup cieutils
class JSONSchemaLoader {
public:
    JSONSchemaLoader();

    explicit JSONSchemaLoader(Ref<const std::filesystem::path> rLibraryRoot);

    explicit JSONSchemaLoader(std::span<const std::filesystem::path> libraryRoots);

    ~JSONSchemaLoader();

    JSONSchemaLoader& operator=(JSONSchemaLoader&& rRHS);

    void addLibraryRoot(Ref<const std::filesystem::path> rLibraryRoot);

    JSONObject load(Ref<const std::string> rURI) const;

private:
    std::vector<std::filesystem::path> _roots;
}; // class JSONSchemaLoader


/// @ingroup cieutils
class JSONSchema {
public:
    JSONSchema();

    JSONSchema(Ref<const JSONObject> rSchema);

    JSONSchema(RightRef<JSONObject> rSchema);

    JSONSchema(
        RightRef<JSONObject> rSchema,
        Ref<const JSONSchemaLoader> rLoader);

    ~JSONSchema();

    JSONSchema& operator=(JSONSchema&& rRHS);

    void resolve(Ref<const JSONSchemaLoader> rLoader = {});

    void validate(Ref<const JSONObject> rJSON) const;

    JSONObject json() const;

    void fillFromDefaults(
        Ref<JSONObject> rJSON,
        Ref<const JSONSchemaLoader> rLoader = {}) const;

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class JSONSchema


/// @ingroup cieutils
std::ostream& operator<<(
    std::ostream& rStream,
    const JSONObject& rJSON);


/// @ingroup cieutils
std::ostream& operator<<(
    Ref<std::ostream> rStream,
    Ref<const JSONSchema> rSchema);


} // namespace cie::io

#include "packages/io/impl/json_impl.hpp"
