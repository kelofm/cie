// --- Utility Includes ---
#include "packages/commandline/inc/KeywordParser.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::utils::detail {


std::string KeywordParser::strip(const std::string& r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::smatch match;
    const bool isAMatch = std::regex_match(
        r_key.begin(),
        r_key.end(),
        match,
        KeywordParser::regex
    );

    CIE_CHECK(
        isAMatch && (match[1].matched != match[2].matched),
        "Invalid key for optional argument: '" + r_key + "'"
    )

    return match[1].matched ? match[1].str() : match[2].str();

    CIE_END_EXCEPTION_TRACING
}


std::string KeywordParser::parse(const std::string& r_key)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return !KeywordParser::strip(r_key).empty() ? r_key : "";

    CIE_END_EXCEPTION_TRACING
}


const std::regex KeywordParser::regex(R"(\-([a-zA-Z_])$|^\-\-([a-zA-Z_][a-zA-Z0-9\-_]*)$)");


} // namespace cie::utils::detail