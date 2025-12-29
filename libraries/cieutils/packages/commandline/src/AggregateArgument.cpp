// --- Utility Includes ---
#include "packages/commandline/inc/AggregateArgument.hpp"
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/commandline/inc/PositionalArgument.hpp"
#include "packages/commandline/inc/KeywordArgument.hpp"
#include "packages/commandline/inc/FlagArgument.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <regex>
#include <sstream>


namespace cie::utils::detail {


class KeyParser
{
public:
    static AbsArgument::Tag parse(const ArgParse::Key& r_key, ArgParse::ArgumentCount nArgs)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        std::match_results<ArgParse::Key::const_iterator> match;
        if (!std::regex_match(r_key.begin(), r_key.end(), match, KeyParser::_regex))
            return AbsArgument::Tag::Invalid;

        if (match[1].matched || match[2].matched)
        {
            if (nArgs == ArgParse::ArgumentCount::None)
                return AbsArgument::Tag::Flag;
            else
                return AbsArgument::Tag::Keyword;
        }
        else if (match[3].matched && nArgs != ArgParse::ArgumentCount::None)
            return AbsArgument::Tag::Positional;

        return AbsArgument::Tag::Invalid;

        CIE_END_EXCEPTION_TRACING
    }

private:
    static const std::regex _regex;
}; // class KeyParser


const std::regex KeyParser::_regex(R"(^(\-[a-zA-Z_])$|^(\-\-[a-zA-Z_][a-zA-Z0-9\-_]*)$|^()$)");


AbsArgument::Tag AbsAggregateArgument::parseTag(const ArgParse::Key& r_key, ArgParse::ArgumentCount nArgs)
{
    return KeyParser::parse(r_key, nArgs);
}


// Check whether the type of the argument is compatible with the expected number of values
void checkArgument(AbsArgument::Tag tag,
                   ArgParse::ArgumentCount nArgs,
                   bool isOptional,
                   const std::string& r_name)
{
    switch (tag)
    {
        case (AbsArgument::Tag::Invalid):
        {
            CIE_THROW(Exception, "Invalid argument '" + r_name + "'")
            break;
        }
        case (AbsArgument::Tag::Positional):
        {
            CIE_CHECK(
                nArgs != ArgParse::ArgumentCount::None,
                "Positional argument '" + r_name + "' must expect at least 1 value"
            )
            CIE_CHECK(
                !isOptional,
                "Positional argument '" + r_name + "' must not be optional"
            )
            break;
        } // Positional
        case (AbsArgument::Tag::Keyword):
        {
            break;
        } // Optional
        case (AbsArgument::Tag::Flag):
        {
            CIE_CHECK(
                nArgs == ArgParse::ArgumentCount::None,
                "Flag '" + r_name + "' must not expect values"
            )
            CIE_CHECK(
                isOptional,
                "Flag '" + r_name + "' must be optional"
            )
            break;
        }
    } // switch tag
}


template <class TArgument>
AggregateArgument<TArgument>::AggregateArgument(ArgumentContainer&& r_arguments,
                                                Size id,
                                                bool isOptional,
                                                ArgParse::DefaultValue&& r_defaultValue,
                                                ArgParse::ArgumentCount nArgs,
                                                const ArgParse::Validator& r_validator,
                                                std::string&& r_docString,
                                                std::string&& r_name)
    : AbsAggregateArgument(),
      IDObject<Size>(id),
      NamedObject(std::move(r_name)),
      _isOptional(isOptional),
      _defaultValue(std::move(r_defaultValue)),
      _nArgs(nArgs),
      _validator(r_validator),
      _docString(std::move(r_docString)),
      _arguments(std::move(r_arguments))
{
    // Validate name
    CIE_CHECK(
        !this->name().empty(),
        "Missing argument name"
    )

    // Check whether the type of the argument is compatible with the expected number of values
    checkArgument(this->tag(), _nArgs, _isOptional, this->name());

    // Validate the number of default values
    if (_nArgs != ArgParse::ArgumentCount::Any && _nArgs != ArgParse::ArgumentCount::NonZero)
    {
        CIE_CHECK(
            r_defaultValue.empty() || int(r_defaultValue.size()) == int(_nArgs),
            (std::stringstream() << "'" << r_name << "' got " << r_defaultValue.size() << " default arguments, but expects " << int(_nArgs)).str()
        )
    }

    // Validate default values
    for (const auto& r_value : r_defaultValue)
        CIE_CHECK(
            _validator(ArgParse::ValueView(&*r_value.begin(), &*r_value.end())),
            "Default value for '" + r_name + "' failed validation"
        )
}


inline void parseOne(const ArgParse::ValueView& r_value,
                     AbsAggregateArgument::Map& r_outputMap,
                     const ArgParse::Validator& r_validator,
                     const std::string& r_argumentName)
{
    ArgParse::Value value(r_value.begin(), r_value.end());

    CIE_CHECK(
        r_validator(r_value),
        "'" + value + "' failed validation for argument '" + r_argumentName + "'"
    )
    r_outputMap.emplace(r_argumentName, std::move(value));
}


template <class TArgument>
bool
AggregateArgument<TArgument>::matchesKey(const ArgParse::KeyView& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (const auto& rArgument : _arguments)
        if (rArgument.matchesKey(r_key))
            return true;

    return false;

    CIE_END_EXCEPTION_TRACING
}


template <class TArgument>
bool
AggregateArgument<TArgument>::matchesKey(const ArgParse::ValueView& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (const auto& rArgument : _arguments)
        if (rArgument.matchesKey(r_key))
            return true;

    return false;

    CIE_END_EXCEPTION_TRACING
}


template <class TArgument>
void
AggregateArgument<TArgument>::parseValues(const ArgParse::ValueViewContainer& r_values, Map& r_outputMap) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    const int nArgs = int(_nArgs);
    if (0 <= nArgs) {
        CIE_CHECK(
            r_values.size() == static_cast<Size>(nArgs),
            "'" << this->name() << "' expected " << int(_nArgs) << " arguments, but got " << r_values.size()
        )

        auto it_value = r_values.begin();
        for (Size index=0; index<static_cast<Size>(nArgs); ++index, ++it_value)
            parseOne(*it_value, r_outputMap, _validator, this->name());
    } else if (_nArgs == ArgParse::ArgumentCount::NonZero) {
        CIE_CHECK(
            0ul < r_values.size(),
            "'" + this->name() + "' expected at least 1 argument, but got " + std::to_string(r_values.size())
        )

        for (const auto& r_value : r_values)
            parseOne(r_value, r_outputMap, _validator, this->name());
    } else if (_nArgs == ArgParse::ArgumentCount::Any) {
        for (const auto& r_value : r_values)
            parseOne(r_value, r_outputMap, _validator, this->name());
    } else
        CIE_THROW(Exception, "Unhandled argument count " + std::to_string(int(_nArgs)))

    CIE_END_EXCEPTION_TRACING
}


template <class TArgument>
const ArgParse::DefaultValue&
AggregateArgument<TArgument>::defaultValue() const
{
    return _defaultValue;
}


template <class TArgument>
bool AggregateArgument<TArgument>::isOptional() const
{
    return _isOptional;
}


template <class TArgument>
ArgParse::ArgumentCount AggregateArgument<TArgument>::nArgs() const
{
    return _nArgs;
}


template <class TArgument>
const std::string& AggregateArgument<TArgument>::docString() const
{
    return _docString;
}


template <class TArgument>
AbsArgument::Tag AggregateArgument<TArgument>::tag() const
{
    return TArgument::tag;
}


template <class TArgument>
const std::string& AggregateArgument<TArgument>::name() const
{
    return NamedObject::name();
}


template <class TArgument>
ArgParse::KeyContainer AggregateArgument<TArgument>::keys() const
{
    ArgParse::KeyContainer keys;
    for (const auto& rArgument : _arguments)
        if (!rArgument.key().empty())
            keys.emplace_back(rArgument.key());
    return keys;
}


// TODO
template <class TStream>
TStream& streamInsert(TStream& rStream, const AbsAggregateArgument& rArgument)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (rArgument.tag() == AbsArgument::Tag::Positional)
        rStream << rArgument.name() << " ";
    else {
        if (rArgument.isOptional())
            rStream << "[";

        for (const auto& r_key : rArgument.keys())
            rStream << r_key << " ";

        if (rArgument.isOptional())
            rStream << "] ";
    }

    const auto nArgs = rArgument.nArgs();
    if (int(nArgs)) {
        rStream << "(";
        if (nArgs == ArgParse::ArgumentCount::Any)
            rStream << "any number";
        else if (nArgs == ArgParse::ArgumentCount::NonZero)
            rStream << "0+";
        else
            rStream << int(nArgs);
        rStream << ")";
    }

    //if (rArgument.tag() != AbsArgument::Tag::Positional)
    //    rStream << rArgument.name() << ", ";

    if (!rArgument.docString().empty())
        rStream << "\n\t" << rArgument.docString();

    if (!rArgument.defaultValue().empty()) {
        rStream << "\n\tDefault: ";
        for (const auto& r_default : rArgument.defaultValue())
            rStream << r_default << " ";
    }

    return rStream;

    CIE_END_EXCEPTION_TRACING
}


OutputStream& operator<<(OutputStream& rStream, const AbsAggregateArgument& rArgument)
{
    return streamInsert(rStream, rArgument);
}


std::ostream& operator<<(std::ostream& rStream, const AbsAggregateArgument& rArgument)
{
    return streamInsert(rStream, rArgument);
}


template class AggregateArgument<PositionalArgument>;


template class AggregateArgument<KeywordArgument>;


template class AggregateArgument<FlagArgument>;


} // namespace cie::utils::detail
