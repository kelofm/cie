#ifndef CIE_UTILS_IO_JSON_HPP
#define CIE_UTILS_IO_JSON_HPP

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
} // namespace cie::concepts::detail



namespace cie::io
{


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
 *  @addtogroup cieutils
 */
class JSONObject
{
public:
    using content_type = nlohmann::json;

private:
    template <class Value>
    class IteratorBase
    {
    public:
        using value_type      = typename CopyConstQualifier<Value,JSONObject>::Type;
        using pointer         = value_type*;
        using reference       = value_type&;
        using difference_type = int;

    public:
        IteratorBase(value_type& r_json);

        IteratorBase(value_type& r_json, difference_type index);

        IteratorBase() = delete;

        IteratorBase(IteratorBase&& r_rhs) = default;

        IteratorBase(const IteratorBase& r_rhs)
        requires concepts::Const<Value> = default;

        IteratorBase& operator=(IteratorBase&& r_rhs) = default;

        IteratorBase& operator=(const IteratorBase& r_rhs)
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

        //bool operator==(const IteratorBase& r_rhs);

        bool operator!=(const IteratorBase& r_rhs);

        bool operator<(const IteratorBase& r_rhs);

        bool operator<=(const IteratorBase& r_rhs);

        bool operator>(const IteratorBase& r_rhs);

        bool operator>=(const IteratorBase& r_rhs);

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

    JSONObject(const JSONObject& r_rhs);

    JSONObject(JSONObject&& r_rhs);

    JSONObject& operator=(const JSONObject& r_rhs);

    JSONObject(const std::string& r_jsonString);

    JSONObject(std::string&& r_jsonString);

    JSONObject(const std::filesystem::path& r_filePath);

    JSONObject(std::istream& r_stream);

    /// Constructor for operator[]
    JSONObject(content_type* p_contents,
               const JSONObject* p_root);

    /// Constructor for operator[] const
    JSONObject(const content_type* p_contents,
               const JSONObject* p_root);

    ~JSONObject();

    /// Get the value associated to the specified key
    JSONObject operator[](const std::string& r_key);

    /// Get the value associated to the specified key
    const JSONObject operator[](const std::string& r_key) const;

    /// Get the specified component if this is an array
    JSONObject operator[](Size index);

    /// Get the specified component if this is an array
    const JSONObject operator[](Size index) const;

    /** Create a new string item with the specified key and value
     *  @details implicitly converts char arrays to std::string
     */
    JSONObject& add(const std::string& r_key,
                    const std::string& r_value,
                    bool allowOverwrite = false);

    /// Create a new item with the specified key and value
    template <concepts::io::SupportedType ValueType>
    JSONObject& add(const std::string& r_key,
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
    bool hasKey(const std::string& r_key) const;

    /// Item size
    Size size() const;

    iterator begin();

    const_iterator begin() const;

    iterator end();

    const_iterator end() const;

    /// Get wrapped object
    const content_type& contents() const;

    /// Get the root that holds the resources
    const JSONObject& root() const;

private:
    /// Set / get helper class
    template <class ValueType>
    struct SetGet
    {
        static void set(JSONObject& r_json,
                         const ValueType& r_value);

        static ValueType as(const JSONObject& r_json);
    };

    /// Type identification and initialization helper
    template <class ValueType>
    struct TypeQuery
    {
        static bool is(const JSONObject& r_json);

        static JSONObject addDefault(JSONObject& r_json,
                                      const std::string& r_key);
    };

    /// Type-safe initialization of a new item's value
    template <concepts::io::SupportedType ValueType>
    JSONObject addDefault(const std::string& r_key);

protected:
    utils::RuntimeConst<content_type> _p_contents;

private:
    const JSONObject* _p_root;
}; // class JSONObject


/// @addtogroup cieutils
std::ostream& operator<<(std::ostream& r_stream, const JSONObject& r_json);


} // namespace cie::io

#include "packages/io/impl/json_impl.hpp"

#endif
