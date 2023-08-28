#ifndef CIE_UTILS_ARG_PARSE_RESULTS_IMPL_HPP
#define CIE_UTILS_ARG_PARSE_RESULTS_IMPL_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/commandline/inc/KeywordParser.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/ArgumentMapper.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/logical.hpp"
#include "packages/compile_time/packages/concepts/inc/functional.hpp"

// --- STL Includes ---
#include <string>
#include <utility>


namespace cie::utils {


namespace detail {
template <class T>
struct ArgParseOptionalTypeConverter
{using Type = T;};

template <>
struct ArgParseOptionalTypeConverter<const char*>
{using Type = std::string;};

template <Size I>
struct ArgParseOptionalTypeConverter<const char (&)[I]>
{using Type = std::string;};

template <>
struct ArgParseOptionalTypeConverter<int>
{using Type = ArgParse::ArgumentCount;};

template <concepts::CallableWith<const ArgParse::ValueView&> TValidator>
struct ArgParseOptionalTypeConverter<TValidator>
{using Type = ArgParse::Validator;};

template <class T, class ...TArgs>
struct ArgParseOptionalTupleConverterImpl
{};

template <Size ...I, class ...TArgs>
struct ArgParseOptionalTupleConverterImpl<std::index_sequence<I...>, TArgs...>
{
private:
    using Original = std::tuple<TArgs...>;

public:
    using Tuple = std::tuple<typename ArgParseOptionalTypeConverter<typename std::tuple_element<I,Original>::type>::Type...>;
};

template <Size S, class ...TArgs>
struct ArgParseOptionalTupleConverter
{using Tuple = typename ArgParseOptionalTupleConverterImpl<std::make_index_sequence<S>, TArgs...>::Tuple;};

template <class ...TArgs>
struct ArgParseOptionalTupleConverter<0,TArgs...>
{using Tuple = std::tuple<>;};
} // namespace detail


inline ArgParse::Results::const_iterator ArgParse::Results::begin() const
{
    return _map.begin();
}


inline ArgParse::Results::const_iterator ArgParse::Results::end() const
{
    return _map.end();
}


inline ArgParse::Results::const_reverse_iterator ArgParse::Results::rbegin() const
{
    return _map.rbegin();
}


inline ArgParse::Results::const_reverse_iterator ArgParse::Results::rend() const
{
    return _map.rend();
}


inline Size ArgParse::Results::size() const
{
    return _map.size();
}


inline bool ArgParse::Results::empty() const
{
    return _map.empty();
}


inline bool ArgParse::Results::has(const Key& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return _map.count(r_key) != 0;

    CIE_END_EXCEPTION_TRACING
}


template <class T>
T ArgParse::Results::get(const Key& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto pair = _map.equal_range(r_key);
    CIE_CHECK(
        pair.first != pair.second,
        *this << std::endl << "No results for key '" << r_key << "'"
    )

    return ArgParse::Results::ValueConverter<T>::convert(pair.first, pair.second);

    CIE_END_EXCEPTION_TRACING
}


template <>
struct ArgParse::Results::ValueConverter<Bool>
{
    static Bool convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        Bool output = false;
        std::stringstream(begin->second) >> std::boolalpha >> output;
        return output;
    }
}; // struct ArgParse::Results::ValueConverter<Bool>


template <>
struct ArgParse::Results::ValueConverter<int>
{
    static int convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return std::stoi(begin->second);
    }
}; // struct ArgParse::Results::ValueConverter<int>


template <>
struct ArgParse::Results::ValueConverter<Size>
{
    static Size convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return std::stoul(begin->second);
    }
}; // struct ArgParse::Results::ValueConverter<Size>


template <>
struct ArgParse::Results::ValueConverter<Float>
{
    static Float convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return std::stof(begin->second);
    }
}; // struct ArgParse::Results::ValueConverter<Float>


template <>
struct ArgParse::Results::ValueConverter<Double>
{
    static Double convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return std::stod(begin->second);
    }
}; // struct ArgParse::Results::ValueConverter<Double>


template <>
struct ArgParse::Results::ValueConverter<std::string>
{
    static std::string convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return begin->second;
    }
}; // struct ArgParse::Results::ValueConverter<std::string>


template <>
struct ArgParse::Results::ValueConverter<std::filesystem::path>
{
    static std::filesystem::path convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return begin->second;
    }
}; // struct ArgParse::Results::ValueConverter<std::filesystem::path>


template <class ...TArgs>
ArgParse& ArgParse::addPositional(std::string&& r_name, TArgs&&... r_optionals)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using OptionalTuple = typename detail::ArgParseOptionalTupleConverter<ct::PackSize<TArgs...>,TArgs...>::Tuple;
    OptionalTuple optionals(std::forward<TArgs>(r_optionals)...);

    auto add = [this, r_name = std::move(r_name)](ArgumentCount argumentCount,
                                                  Validator r_validator,
                                                  DefaultValue r_defaultValue,
                                                  std::string r_docString) mutable -> ArgParse&
    {
        return this->addArgument(std::move(r_name),
                                 {""},
                                 false,
                                 argumentCount,
                                 r_validator,
                                 std::move(r_defaultValue),
                                 std::move(r_docString));
    };

    return ct::ArgumentMapper<ArgumentCount,Validator,DefaultValue,std::string>(
        ArgumentCount::One,
        defaultValidator,
        DefaultValue {},
        std::string("")
    ).mapTuple<ArgParse&>(add, std::move(optionals));

    CIE_END_EXCEPTION_TRACING
}


template <class ...TArgs>
ArgParse& ArgParse::addKeyword(KeyContainer&& r_keys, TArgs&&... r_optionals)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using OptionalTuple = typename detail::ArgParseOptionalTupleConverter<ct::PackSize<TArgs...>,TArgs...>::Tuple;
    OptionalTuple optionals(std::forward<TArgs>(r_optionals)...);

    auto add = [this, r_keys = std::move(r_keys)](bool isOptional,
                                                  ArgumentCount argumentCount,
                                                  Validator r_validator,
                                                  DefaultValue r_defaultValue,
                                                  std::string r_docString) mutable -> ArgParse&
    {
        std::string name;
        if (!r_keys.empty())
            name = detail::KeywordParser::strip(r_keys.front());
        else
            CIE_THROW(Exception, "Attempting to add a keyword argument without keys assigned to it.")

        return this->addArgument(std::move(name),
                                 std::move(r_keys),
                                 isOptional,
                                 argumentCount,
                                 r_validator,
                                 std::move(r_defaultValue),
                                 std::move(r_docString));
    };

    return ct::ArgumentMapper<bool,ArgumentCount,Validator,DefaultValue,std::string>(
        true,
        ArgumentCount::One,
        defaultValidator,
        DefaultValue {},
        std::string("")
    ).mapTuple<ArgParse&>(add, std::move(optionals));

    CIE_END_EXCEPTION_TRACING
}


template <class ...TArgs>
ArgParse& ArgParse::addFlag(KeyContainer&& r_keys, TArgs&&... r_optionals)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using OptionalTuple = typename detail::ArgParseOptionalTupleConverter<ct::PackSize<TArgs...>,TArgs...>::Tuple;
    OptionalTuple optionals(std::forward<TArgs>(r_optionals)...);

    auto add = [this, r_keys = std::move(r_keys)](std::string r_docString) mutable -> ArgParse&
    {
        std::string name;
        if (!r_keys.empty())
            name = detail::KeywordParser::strip(r_keys.front());
        else
            CIE_THROW(Exception, "Attempting to add a keyword argument without keys assigned to it.")

        return this->addArgument(std::move(name),
                                 std::move(r_keys),
                                 true,
                                 ArgumentCount::None,
                                 this->defaultValidator,
                                 {},
                                 std::move(r_docString));
    };

    return ct::ArgumentMapper<std::string>(std::string("")).mapTuple<ArgParse&>(add, std::move(optionals));

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils


#endif