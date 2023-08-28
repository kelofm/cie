// --- Utility Includes ---
#include "packages/commandline/inc/FlagArgument.hpp"
#include "packages/commandline/inc/KeywordParser.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils::detail {


FlagArgument::FlagArgument(const ArgParse::Key& r_key)
    : AbsArgument(r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::smatch match;
    CIE_CHECK(
        std::regex_match(r_key.begin(), r_key.end(), match, KeywordParser::regex),
        "Invalid key: '" + r_key + "'"
    )

    if (match[1].matched) // single character flag
        _regex = std::regex(R"(^\-[a-zA-Z]*)" + match[1].str() + R"([a-zA-Z]*$)");
    else if (match[2].matched)
        _regex = std::regex(R"(^\-\-)" + match[2].str() + R"($)");

    CIE_END_EXCEPTION_TRACING
}


bool FlagArgument::matchesKey(const ArgParse::KeyView& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::match_results<ArgParse::KeyView::const_iterator> match;
    return std::regex_match(
        r_key.begin(),
        r_key.end(),
        match,
        _regex
    );

    CIE_END_EXCEPTION_TRACING
}


bool FlagArgument::matchesKey(const ArgParse::ValueView& r_key) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::match_results<decltype(r_key.begin())> match;
    return std::regex_match(
        r_key.begin(),
        r_key.end(),
        match,
        _regex
    );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils::detail