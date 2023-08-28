#ifndef CIE_UTILS_COMMAND_LINE_KEYWORD_PARSER_HPP
#define CIE_UTILS_COMMAND_LINE_KEYWORD_PARSER_HPP

// --- STL Includes ---
#include <regex>


namespace cie::utils::detail {


struct KeywordParser
{
public:
    static std::string strip(const std::string& r_key);

    static std::string parse(const std::string& r_key);

    static const std::regex regex;
};


} // namespace cie::utils::detail


#endif