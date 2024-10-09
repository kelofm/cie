#ifndef CIE_UTILS_ARG_PARSE_HPP
#define CIE_UTILS_ARG_PARSE_HPP

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/logging/inc/OutputStream.hpp"
#include "packages/types/inc/NamedObject.hpp"

// --- STL Includes ---
#include <functional>
#include <ranges>
#include <string_view>
#include <map>
#include <ostream>


namespace cie::utils {


namespace detail {
// Required for friend declaration in a different namespace
class AbsAggregateArgument;
} // namespace detail


///@ingroup cieutils
class ArgParse : utils::NamedObject
{
public:
    CIE_DEFINE_CLASS_POINTERS(ArgParse)

    using Key = std::string;

    using KeyView = std::basic_string_view<Key::value_type>;

    using KeyContainer = DynamicArray<Key>;

    using Value = std::string;

    using ValueContainer = DynamicArray<Value>;

    using ValueView = std::ranges::subrange<const Value::value_type*>;

    using ValueViewContainer = DynamicArray<ValueView>;

    using DefaultValue = ValueContainer;

    using Validator = std::function<bool(const ValueView&)>;

    class ArgumentCount
    {
    public:
        static const ArgumentCount None;

        static const ArgumentCount One;

        static const ArgumentCount NonZero;

        static const ArgumentCount Any;

        static const ArgumentCount exactly(Size count);

    public:
        constexpr ArgumentCount()
            : ArgumentCount(1)
        {}

        constexpr ArgumentCount(int value)
            : _value(value)
        {}

        constexpr ArgumentCount(ArgumentCount&& r_rhs) = default;

        constexpr ArgumentCount(const ArgumentCount& r_rhs) = default;

        constexpr ArgumentCount& operator=(ArgumentCount&& r_rhs) = default;

        constexpr ArgumentCount& operator=(const ArgumentCount& r_rhs) = default;

        friend bool operator==(const ArgumentCount lhs, const ArgumentCount rhs)
        {return lhs._value == rhs._value;}

        friend bool operator==(const ArgumentCount lhs, const int rhs)
        {return lhs._value == rhs;}

        friend bool operator==(const int lhs, const ArgumentCount rhs)
        {return lhs == rhs._value;}

        friend bool operator!=(const ArgumentCount lhs, const ArgumentCount rhs)
        {return lhs._value != rhs._value;}

        friend bool operator!=(const ArgumentCount lhs, const int rhs)
        {return lhs._value != rhs;}

        friend bool operator!=(const int lhs, const ArgumentCount rhs)
        {return lhs != rhs._value;}

        explicit operator int () const
        {return _value;}

    private:
        int _value;
    }; // class ArgumentCount

    class Results
    {
    public:
        using Key = ArgParse::Key;

        CIE_DEFINE_CLASS_POINTERS(Results)

        using Value = ArgParse::Value;

    private:
        using Map = std::multimap<Key,Value>;

        friend class detail::AbsAggregateArgument;

        friend class ArgParse;

    public:
        using const_iterator = Map::const_iterator;

        using const_reverse_iterator = Map::const_reverse_iterator;

        using ValueIterator = const_iterator;

    public:
        CIE_DEFINE_CLASS_DEFAULTS(Results)

        const_iterator begin() const;

        const_iterator end() const;

        const_reverse_iterator rbegin() const;

        const_reverse_iterator rend() const;

        Size size() const;

        bool empty() const;

        /// Check whether the internal container has a value for the specified key.
        bool has(const Key& r_key) const;

        std::ranges::subrange<const_iterator> at(const Key& r_key) const;

        std::ranges::subrange<const_iterator> operator[](const Key& r_key) const;

        /** Search for a key and convert its associated value to the specified type.
        *
        *  @param r_key: key to be searched for in the internal container. An exception is thrown
        *                if the key is not found.
        *  @details This function redirectrects to @ref ArgParse::Results::ValueConverter::convert,
        *  which has specializations for a number of types. Missing specializations
        *  for the specified type will result in a linker error. Additional specializations
        *  can be provided in this case.
        *  For implemented specializations, see @ref ArgParse::Results::ValueConverter::convert.
        */
        template <class T>
        T get(const Key& r_key) const;

        /** Convert a value to the desired type.
        *
        *  Overload this class to add conversion rules to a specific type.
        *  Implemented overloads:
        *  - @ref Bool
        *  - int
        *  - @ref Size
        *  - @ref Float
        *  - @ref Double
        *  - std::string
        */
        template <class T>
        struct ValueConverter
        {
            static T convert(ValueIterator begin, ValueIterator end);
        }; // struct ValueConverter

        OutputStream& serialize(OutputStream& r_stream) const;

    private:
        Results(Map&& r_map);

        Results(const Map& r_map);

    private:
        Map _map;
    }; // class Results

    /// Always true.
    static const Validator defaultValidator;

    /// Checks whether the input string is a path that points to an existing file or directory.
    static const Validator pathValidator;

    /// Checks whether the input string is a path that points to a directory.
    static const Validator directoryPathValidator;

    // Checks whether the input string is a path that points to a file.
    static const Validator filePathValidator;

public:
    /// @brief Default constructor initializing pimpl.
    ArgParse();

    ArgParse(std::string&& r_name);

    ArgParse(const std::string& r_name);

    ArgParse(ArgParse&& r_rhs) = default;

    ArgParse& operator=(ArgParse&& r_rhs) = default;

    /// @brief Explicit constructor required by pimpl.
    ~ArgParse();

    /// @brief Check whether an argument with the specified name exists
    bool has(const Key& r_key) const;

    /** @brief Add a positional argument.
     *
     *  @param r_name name of the argument. Must not be empty.
     *  @param r_optionals optional arguments (in any order):
     *                     - @p nArgs number of expected values. Must not be @ref ArgumentCount::None but may
     *                                be variable (@ref ArgumentCount::Any or @ref ArgumentCount::NonZero).
     *                     - @p r_validator functor that validates each parsed value.
     *                     - @p r_docString description/help for the variable.
     *
     *  @details A positional argument is a required argument that expects at least one value before
     *           any keywords or flags are specified. Since they are not optional, positional arguments
     *           have no default arguments.
     */
    template <class ...TArgs>
    ArgParse& addPositional(std::string&& r_name, TArgs&&... r_optionals);

    /** @brief Add a keyword argument.
     *
     *  @param r_keys container of keys with at least one item.
     *  @param r_optionals optional arguments (in any order):
     *                     - @p isOptional indicates whether the argument is optional.
                                        If false, @ref parseArguments throws an error if the argument
     *                                     is not found in the input.
     *                     - @p nArgs indicates how many values the argument expects. Pass @ref ArgumentCount::Any
     *                                or @ref ArgumentCount::NonZero for a variable number of values.
     *                     - @p r_validator functor that validates each parsed value.
     *                     - @p r_defaultValue container of default values which are used if the argument
     *                                         is not found during parsing. Non-optional arguments must not
     *                                         have default values. If the argument is optional, the number
     *                                         of default values must match @p nArgs.
     *                     - @p r_docString description/help for the variable.
     *
     *  @details A keyword is an optional or required argument that is identified by at least one key,
     *           and may expect any number of values except 0. Values of keyword arguments must follow
     *           their respective keys. All passed keys must match one of two key patterns:
     *           - beginning with '-' followed by a single English character [a-zA-Z]
     *           - beginning with '--' followed by a single English character [a-zA-Z] then any number of alphanumeric characters [a-zA-Z0-9]
     */
    template <class ...TArgs>
    ArgParse& addKeyword(KeyContainer&& r_keys, TArgs&&... r_optionals);

    /** @brief Add a flag.
     *
     *  @param r_keys a container of keys with at least one item.
     *  @param r_optionals optional arguments (in any order):
     *                     - @p r_docString description/help for the variable.
     *
     *  @details A flag is an optional variable that must not expect explicit values,
     *           but has a binary true/false value depending whether it is set or not respectively.
     *           It must have at least one key, and all passed keys must match one of two key patterns:
     *           - beginning with '-' followed by a single English character [a-zA-Z]
     *           - beginning with '--' followed by a single English character [a-zA-Z] then any number of alphanumeric characters [a-zA-Z0-9]
     */
    template <class ...TArgs>
    ArgParse& addFlag(KeyContainer&& r_keys, TArgs&&... r_optionals);

    /** @brief Deduce the argument type and add it to the list
     *
     *  @details Positional argument:
     *           - has exactly one key, and that is empty
     *           - is not optional
     *           - has no default values
     *           - can have any number of values excefpt @ref ArgumentCount::None
     *
     *           Keyword Argument:
     *               - has at least one key, and it must match a key pattern
     *               - can be optional
     *               - may have default values
     *               - can have any number of values except @ref ArgumentCount::None
     *
     *           Flag argument:
     *               - has at least one key, and it must match a key pattern
     *               - can be optional
     *               - does not have a default value (but has a built-in false value if not set)
     *               - does not have any arguments
     *
     *  @param r_name unique name of the argument
     *  @param r_keys container of keys
     *  @param isOptional indicates whether the argument is optional.
     *                    If false, @ref parseArguments throws an error if the argument
     *                    is not found in the input.
     *  @param nArgs indicates how many values the argument expects. Pass @ref ArgumentCount::Any
     *               or @ref ArgumentCount::NonZero for a variable number of values, in which case
     *               values will be consumed greedily during parsing. Only one positional argument
     *               is allowed to have a variable number of arguments.
     *  @param r_validator functor that validates each parsed value
     *  @param r_defaultValue container of default values which are used if the argument
     *                        is not found during parsing. Non-optional arguments must not
     *                        have default values. If the argument is optional, the number
     *                        of default values must match @p nArgs.
     *  @param r_docString description/help for the argument
     */
    ArgParse& addArgument(std::string&& r_name,
                          KeyContainer&& r_keys,
                          bool isOptional,
                          ArgumentCount nArgs,
                          const Validator& r_validator,
                          DefaultValue&& r_defaultValue,
                          std::string&& r_docString);

    /** @brief Parse the input arguments and return a map containing all parsed and default values
     *
     *  @param argc number of arguments
     *  @param argv array of arguments of size 'argc'
     *  @return { @p argument_name, @p argument_values }
     *
     *  @note The returned map contains every registered argument.
     */
    Results parseArguments(int argc, const char* argv[]) const;

    /// @brief Push help string to the provided stream.
    OutputStream& help(OutputStream& r_stream) const;

    /// @brief Push help string to the provided stream.
    std::ostream& help(std::ostream& r_stream) const;

private:
    ArgParse(const ArgParse& r_rhs) = delete;

    ArgParse& operator=(const ArgParse& r_rhs) = delete;

public:
    /** @brief Return a regex-based validator checking basic types.
     *
     * @details Explicit specializations are implemented for the following classes/built-ins:
     *          - int
     *          - @ref Size
     *          - @ref Double
     */
    template <class T>
    static Validator validatorFactory();

private:
    class Impl;
    std::unique_ptr<Impl> _p_impl;
}; // class ArgParse


OutputStream& operator<<(OutputStream& r_stream, const ArgParse::Results& r_parsedArguments);


std::ostream& operator<<(std::ostream& r_stream, const ArgParse::Results& r_parsedArguments);


} // namespace cie::utils

#include "packages/commandline/impl/ArgParse_impl.hpp"

#endif
