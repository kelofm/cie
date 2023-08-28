#ifndef CIE_UTILS_COMMAND_LINE_FILE_SYSTEM_EXTENSION_IMPL_HPP
#define CIE_UTILS_COMMAND_LINE_FILE_SYSTEM_EXTENSION_IMPL_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::utils {


template <>
struct ArgParse::Results::ValueConverter<std::filesystem::path>
{
    static std::filesystem::path convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        CIE_CHECK(std::distance(begin, end) == 1, "Expecting exactly 1 string, but got " << std::distance(begin, end))
        return std::filesystem::path(begin->second.begin(), begin->second.end());
    }
};


} // namespace cie::utils


#endif