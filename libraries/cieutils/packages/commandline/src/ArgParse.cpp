// --- Utility Incldues ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/commandline/inc/AbsArgument.hpp"
#include "packages/commandline/inc/AggregateArgument.hpp"
#include "packages/commandline/inc/PositionalArgument.hpp"
#include "packages/commandline/inc/KeywordArgument.hpp"
#include "packages/commandline/inc/FlagArgument.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/logging/inc/StreamInterceptor.hpp"

// --- STL Includes ---
#include <unordered_set>
#include <algorithm>
#include <ranges>
#include <sstream>
#include <limits>
#include <iostream>
#include <regex>
#include <filesystem>


namespace cie::utils {


const ArgParse::ArgumentCount ArgParse::ArgumentCount::None(0);


const ArgParse::ArgumentCount ArgParse::ArgumentCount::One(1);


const ArgParse::ArgumentCount ArgParse::ArgumentCount::NonZero(-1);


const ArgParse::ArgumentCount ArgParse::ArgumentCount::Any(-2);


std::ranges::subrange<ArgParse::Results::const_iterator> ArgParse::Results::at(const Key& r_key) const
{
    return this->operator[](r_key);
}


std::ranges::subrange<ArgParse::Results::const_iterator> ArgParse::Results::operator[](const Key& r_key) const
{
    const auto pair = _map.equal_range(r_key);
    return {pair.first, pair.second};
}


ArgParse::Results::Results(Map&& r_map)
    : _map(std::move(r_map))
{
}


ArgParse::Results::Results(const Map& r_map)
    : _map(r_map)
{
}


const ArgParse::Validator ArgParse::defaultValidator = [](const auto&){return true;};


const ArgParse::Validator ArgParse::pathValidator = [](const auto& r_value)
{return std::filesystem::exists(std::filesystem::path(r_value.begin(), r_value.end()));};


const ArgParse::Validator ArgParse::directoryPathValidator = [](const auto& r_value)
{return std::filesystem::is_directory(std::filesystem::path(r_value.begin(), r_value.end()));};


const ArgParse::Validator ArgParse::filePathValidator = [](const auto& r_value)
{return std::filesystem::is_regular_file(std::filesystem::path(r_value.begin(), r_value.end()));};


class ArgParse::Impl
{
public:
    using Argument = detail::AbsAggregateArgument;

    using ArgumentContainer = DynamicArray<Argument::SharedPointer>;

    ~Impl();

    ArgumentContainer _arguments;
}; // class ArgParse::Impl


ArgParse::Impl::~Impl()
{
}


ArgParse::ArgParse()
    : utils::NamedObject(""),
      _p_impl(new ArgParse::Impl)
{
}


ArgParse::ArgParse(std::string&& r_name)
    : utils::NamedObject(std::move(r_name)),
      _p_impl(new ArgParse::Impl)
{
}


ArgParse::ArgParse(const std::string& r_name)
    : utils::NamedObject(r_name),
      _p_impl(new ArgParse::Impl)
{
}


ArgParse::~ArgParse()
{
}


bool ArgParse::has(const Key& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (const auto& rp_argument : this->_p_impl->_arguments)
        if (rp_argument->matchesKey(r_key))
            return true;
    return false;

    CIE_END_EXCEPTION_TRACING
}


ArgParse& ArgParse::addArgument(std::string&& r_name,
                                KeyContainer&& r_keys,
                                bool isOptional,
                                ArgumentCount nArgs,
                                const Validator& r_validator,
                                DefaultValue&& r_defaultValue,
                                std::string&& r_docString)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check keys
    CIE_CHECK(
        !r_keys.empty(),
        "At least one key must be provided for each argument. Positional arguments are represented by empty strings for keys."
    )

    // Check argument name
    CIE_CHECK(
        !r_name.empty(),
        "Empty argument name"
    )

    for (const auto& rp_argument : _p_impl->_arguments)
    {
        CIE_CHECK(
            rp_argument->name() != r_name,
            "Argument with name '" + r_name + "' already exists"
        )
    }

    // Check whether the new argument's keys clash with existing ones
    const auto& r_arguments = _p_impl->_arguments;
    for (const auto& rp_argument : r_arguments)
        for (const auto& r_key : r_keys)
    CIE_CHECK(
        !rp_argument->matchesKey(r_key),
        "Key '" + r_key + "' of new argument '" + r_name + "' clashes with argument '" + rp_argument->name() + "'"
    )

    // Init
    const Size id = _p_impl->_arguments.size();
    auto tag = detail::AbsArgument::Tag::Invalid;

    // Get tag
    for (const auto& r_key : r_keys)
    {
        const auto localTag = Impl::Argument::parseTag(r_key, nArgs);
        CIE_CHECK(
            localTag != detail::AbsArgument::Tag::Invalid,
            "Invalid key '" + r_key + "'"
        )

        if (tag == detail::AbsArgument::Tag::Invalid)
            tag = localTag;
        else if (tag != localTag)
            CIE_THROW(Exception, "Argument type mismatch for key'" + r_key + "' in argument '" + r_name + "'")
    }

    // Construct arguments
    switch (tag)
    {
        case (detail::AbsArgument::Tag::Invalid):
        {
            CIE_THROW(Exception, "Invalid argument '" + r_name + "'")
            break;
        } // Invlaid
        case (detail::AbsArgument::Tag::Positional):
        {
            CIE_CHECK(
                r_keys.size() == 1,
                "Positional argument '" + r_name + "' must have exactly 1 key (and it must be empty)"
            )
            const auto& r_key = r_keys.front();
            CIE_CHECK(
                r_key.empty(),
                "Positional argument '" + r_name + "' has a non-empty key"
            )
            CIE_CHECK(
                r_defaultValue.empty(),
                "Positional argument '" + r_name + "' got default values"
            )

            // Check positional arguments' compatibility
            for (const auto& rp_argument : std::ranges::subrange(_p_impl->_arguments.rbegin(), _p_impl->_arguments.rend()))
                if (rp_argument->tag() == detail::AbsArgument::Tag::Positional)
                {
                    CIE_CHECK(
                        rp_argument->nArgs() != ArgumentCount::Any && rp_argument->nArgs() != ArgumentCount::NonZero,
                        "Attempt to add positional argument '" + r_name + "' after another positional argument '" + rp_argument->name() + "' with variable number of values"
                    )
                    break;
                }

            using ArgumentType = detail::PositionalArgument;
            using AggregateType = detail::AggregateArgument<ArgumentType>;
            AggregateType::ArgumentContainer arguments;
            arguments.emplace_back();

            _p_impl->_arguments.emplace_back(new AggregateType(
                std::move(arguments),
                id,
                isOptional,
                std::move(r_defaultValue),
                nArgs,
                r_validator,
                std::move(r_docString),
                std::move(r_name)
            ));
            break;
        } // Positional
        case (detail::AbsArgument::Tag::Keyword):
        {
            using ArgumentType = detail::KeywordArgument;
            using AggregateType = detail::AggregateArgument<ArgumentType>;
            AggregateType::ArgumentContainer arguments;
            arguments.reserve(r_keys.size());

            for (const auto& r_key : r_keys)
                arguments.emplace_back(r_key);

            _p_impl->_arguments.emplace_back(new AggregateType(
                std::move(arguments),
                id,
                isOptional,
                std::move(r_defaultValue),
                nArgs,
                r_validator,
                std::move(r_docString),
                std::move(r_name)
            ));
            break;
        } // Keyword
        case (detail::AbsArgument::Tag::Flag):
        {
            using ArgumentType = detail::FlagArgument;
            using AggregateType = detail::AggregateArgument<ArgumentType>;
            AggregateType::ArgumentContainer arguments;
            arguments.reserve(r_keys.size());

            for (const auto& r_key : r_keys)
                arguments.emplace_back(r_key);

            _p_impl->_arguments.emplace_back(new AggregateType(
                std::move(arguments),
                id,
                isOptional,
                std::move(r_defaultValue),
                nArgs,
                r_validator,
                std::move(r_docString),
                std::move(r_name)
            ));
            break;
        } // Flag
    } // switch tag

    return *this;

    CIE_END_EXCEPTION_TRACING
}


ArgParse::ValueView CStringToView(const char* begin)
{
    CIE_BEGIN_EXCEPTION_TRACING

    const char* end = begin;
    while (*end != '\0')
        ++end;

    return ArgParse::ValueView(begin, end);

    CIE_END_EXCEPTION_TRACING
}


// Collect values from the input until the required amount is reached
// or, if the argument expects a variable number of values, until a
// value matches one of the arguments. In the latter case, return the
// a pointer to the matched variable, otherwise return nullptr.
detail::AbsAggregateArgument::SharedConstPointer
consumeValues(int& argc,
              const char**& argv,
              const DynamicArray<detail::AbsAggregateArgument::SharedPointer> r_arguments,
              ArgParse::ValueViewContainer& r_views,
              int maxValues)
{
    CIE_BEGIN_EXCEPTION_TRACING

    detail::AbsAggregateArgument::SharedConstPointer p_matched;

    // Loop through values until:
    //  - end of values
    //  - consumed a fixed number of values (maxValues)
    //  - a value matches at least one of the arguments
    for (; 0 < argc && 0 < maxValues && !p_matched; --argc, --maxValues, ++argv)
    {
        // Create a view from the current value
        auto view = CStringToView(*argv);

        // If this value matches with any of the arguments,
        // it means that we've run out of values, or the input
        // is invalid (this case would handled caught later on).
        const auto it_argument = std::find_if(
            r_arguments.begin(),
            r_arguments.end(),
            [&view](const auto& rp_arg){return rp_arg->matchesKey(view);}
        );

        if (it_argument != r_arguments.end())
        {
            p_matched = *it_argument;
            break;
        }

        // Add a view of the current value
        r_views.emplace_back(std::move(view));
    }

    return p_matched;

    CIE_END_EXCEPTION_TRACING
}


ArgParse::Results ArgParse::parseArguments(int argc, const char* argv[]) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    ArgParse::Results::Map outputMap;

    // Create a map that tracks whether an argument
    // appeared in the input
    std::unordered_map<Impl::Argument::SharedConstPointer,bool> argumentFoundMap;
    for (const auto& rp_argument : _p_impl->_arguments)
        argumentFoundMap.emplace(rp_argument, false);

    // Collect positional arguments
    DynamicArray<Impl::Argument::SharedConstPointer> positionals;
    {
        int expectedPositionalCount = 0;

        for (const auto& rp_argument : _p_impl->_arguments)
            if (rp_argument->tag() == detail::AbsArgument::Tag::Positional)
            {
                positionals.emplace_back(rp_argument);
                const auto nArgs = rp_argument->nArgs();
                if (nArgs == ArgumentCount::Any || nArgs == ArgumentCount::NonZero)
                {
                    expectedPositionalCount = argc;
                    break;
                } // variable argument count
                else
                {
                    expectedPositionalCount += int(nArgs);
                } // Exact argument count
            } // if positional argument

        // Check whether the expected number of positional values is valid
        CIE_CHECK(
            expectedPositionalCount <= argc,
            "Expecting " << expectedPositionalCount << " positional argument" << (1 < expectedPositionalCount ? "s" : "") << ", but " << argc << " were given"
        )
    } // construct separators

    // Assign values to positional arguments
    ValueViewContainer views;

    for (const auto& rp_positional : positionals)
    {
        const auto nArgs = rp_positional->nArgs();
        views.clear();

        // Consume values
        // Note: argc and argv are updated as values are consumed
        auto p_matched = consumeValues(
            argc,
            argv,
            _p_impl->_arguments,
            views,
            nArgs != ArgumentCount::Any && nArgs != ArgumentCount::NonZero ? int(nArgs) : argc
        );

        // A value matched an argument => end of positionals or error
        if (p_matched)
            break;

        // Pass the collected views to the positional argument
        // which will validate them. Then push the validated
        // values to the output map.
        try {rp_positional->parseValues(views, outputMap);}
        catch (Exception& r_exception)
        {
            std::stringstream stream;
            this->help(stream);
            CIE_RETHROW(r_exception, stream.str())
        }
        argumentFoundMap[rp_positional] = true;
    }

    // Assign values to keyword and flag arguments
    while (0 < argc)
    {
        // Assume the current value is a key

        // Get argument matching the key.
        // Possible results:
        //  - 'views' is empty and 'p_matched' is valid:
        //    argument is found and no additional values are consumed (OK)
        //
        //  - 'views' is not empty and 'p_matched' is valid:
        //    unexpected positional arguments were found before finding a keyword (INPUT ERROR)
        //
        //  - 'views' is empty and 'p_matched' is invalid:
        //    this should not happen and means there is a bug somewhere (ERROR)
        //
        //  - 'views' is not empty and 'p_matched' is invalid
        //    unexpected trailing values (INPUT ERROR)
        views.clear();
        auto p_matched = consumeValues(
            argc,
            argv,
            _p_impl->_arguments,
            views,
            argc
        );

        if (!views.empty())
        {
            std::stringstream stream;
            stream << "Unexpected arguments or unrecognized keys: ";
            for (const auto& r_view : views)
                stream << "'" << std::string(r_view.begin(), r_view.end()) << "' ";
            CIE_THROW(Exception, stream.str())
        }

        CIE_CHECK(
            p_matched,
            "Expecting a key but found none"
        )

        // From here on, 'p_matched' is valid and 'views' is empty.
        // => move on to the values (if any)
        --argc;
        ++argv;

        // Collect values
        const auto nArgs = p_matched->nArgs();
        consumeValues(
            argc,
            argv,
            _p_impl->_arguments,
            views,
            nArgs != ArgumentCount::Any && nArgs != ArgumentCount::NonZero ? int(nArgs) : argc
        );

        // Check name collision
        if (outputMap.contains(p_matched->name()))
            CIE_THROW(Exception, "Duplicate key: '" + p_matched->name() + "'")

        // Pass the collected views to the argument,
        // which validates and registers them.
        try {p_matched->parseValues(views, outputMap);}
        catch (Exception& r_exception)
        {
            std::stringstream stream;
            this->help(stream);
            CIE_RETHROW(r_exception, stream.str())
        }

        argumentFoundMap[p_matched] = true;
    } // for i_argument < argc

    // Check whether all required arguments were found
    // and assign default values to arguments that were not found.
    for (const auto& r_pair : argumentFoundMap)
    {
        if (!r_pair.second)
        {
            if (!r_pair.first->isOptional())
                CIE_THROW(Exception, "Required argument '" + r_pair.first->name() +"' was not found")
            else if (r_pair.first->tag() == detail::AbsArgument::Tag::Flag)
                outputMap.emplace(r_pair.first->name(), "false");
            else
                for (const auto& r_defaultValue : r_pair.first->defaultValue())
                    outputMap.emplace(r_pair.first->name(), r_defaultValue);
        }
        else if (r_pair.first->tag() == detail::AbsArgument::Tag::Flag)
            outputMap.emplace(r_pair.first->name(), "true");
    }

    return ArgParse::Results(std::move(outputMap));

    CIE_END_EXCEPTION_TRACING
}


template <class TStream, class TArguments>
TStream& helpImpl(TStream& r_stream, const std::string& r_header, const TArguments& r_arguments)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_stream << r_header;
    for (const auto& rp_argument : r_arguments)
        r_stream << *rp_argument << std::endl;

    return r_stream;

    CIE_END_EXCEPTION_TRACING
}


OutputStream& ArgParse::help(OutputStream& r_stream) const
{
    return helpImpl(r_stream,
                    this->name().empty() ? "" : (this->name() + ":\n"),
                    _p_impl->_arguments);
}


std::ostream& ArgParse::help(std::ostream& r_stream) const
{
    return helpImpl(r_stream,
                    this->name().empty() ? "" : (this->name() + ":\n"),
                    _p_impl->_arguments);
}


struct RegexBasedValidator
{
    inline static bool isValid(const ArgParse::ValueView& r_value, const std::regex& r_regex)
    {
        CIE_BEGIN_EXCEPTION_TRACING
        std::match_results<ArgParse::KeyView::const_iterator> match;
        return std::regex_match(r_value.begin(), r_value.end(), match, r_regex);
        CIE_END_EXCEPTION_TRACING
    }
}; // class RegexBasedValidator


const std::regex integerRegex(R"(^0|(?:-?[1-9]+[0-9]*)$)");


const std::regex unsignedIntegerRegex(R"(^0|(?:[1-9]+[0-9]*)$)");


const std::regex floatingPointRegex(R"(^-?(?:(?:(?:[1-9][0-9]*)(?:\.[0-9]*)?)|(?:0(?:\.[0-9]*)?))(?:[eE][\+-]?[0-9]+)?$)");


template <>
ArgParse::Validator ArgParse::validatorFactory<int>()
{
    return [](const auto& r_value) {return RegexBasedValidator::isValid(r_value, integerRegex);};
}


template <>
ArgParse::Validator ArgParse::validatorFactory<Size>()
{
    return [](const auto& r_value) {return RegexBasedValidator::isValid(r_value, unsignedIntegerRegex);};
}


template <>
ArgParse::Validator ArgParse::validatorFactory<Double>()
{
    return [](const auto& r_value) {return RegexBasedValidator::isValid(r_value, floatingPointRegex);};
}


OutputStream& ArgParse::Results::serialize(OutputStream& r_stream) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    constexpr const auto indent = "    ";

    r_stream << "{";

    if (!this->empty()) [[likely]]
        r_stream << std::endl;

    ArgParse::Key key = "";
    for (const auto& r_pair : *this)
    {
        if (key == r_pair.first)
            r_stream << ", ";
        else
        {
            if (!key.empty()) [[likely]]
                r_stream << "}" << std::endl;

            key = r_pair.first;
            r_stream << indent << key << " : {";
        }

        r_stream << r_pair.second;
    }

    if (!this->empty()) [[likely]]
        r_stream << "}" << std::endl;

    r_stream << "}";

    return r_stream;

    CIE_END_EXCEPTION_TRACING
}


OutputStream& operator<<(OutputStream& r_stream, const ArgParse::Results& r_parsedResults)
{
    return r_parsedResults.serialize(r_stream);
}


std::ostream& operator<<(std::ostream& r_stream, const ArgParse::Results& r_parsedResults)
{
    CIE_BEGIN_EXCEPTION_TRACING

    OutputStream::SharedPointer p_wrapper(new OutputStream(&r_stream));
    std::stringstream stringStream;
    OutputStream::SharedPointer p_stream(new OutputStream(&stringStream));
    StreamInterceptor redirector(p_stream, p_wrapper);
    r_parsedResults.serialize(*p_stream);

    return r_stream;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils