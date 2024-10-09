#ifndef CIE_UTILS_COMMAND_LINE_FILE_SYSTEM_EXTENSION_HPP
#define CIE_UTILS_COMMAND_LINE_FILE_SYSTEM_EXTENSION_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::utils {


/// @ingroup cieutils
template <>
struct ArgParse::Results::ValueConverter<std::filesystem::path>;


} // namespace cie::utils

#include "packages/commandline/impl/file_system_extension_impl.hpp"

#endif
