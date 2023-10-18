#ifndef CIE_UTILS_LOGGING_LOGGER_SINGLETON_HPP
#define CIE_UTILS_LOGGING_LOGGER_SINGLETON_HPP

// --- Internal Includes ---
#include "packages/logging/inc/Logger.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::utils {


namespace detail {
extern const std::filesystem::path DEFAULT_LOGGER_FILE_PATH;
} // namespace detail


/// @addtogroup cieutils
class LoggerSingleton
{
public:
    static Logger& get();

    static Logger& get(const std::filesystem::path& r_filepath);

    static LoggerPtr getPtr();

    static LoggerPtr getPtr(const std::filesystem::path& r_filePath);

private:
    static void initialize(const std::filesystem::path& r_filePath = detail::DEFAULT_LOGGER_FILE_PATH);

private:
    static LoggerPtr _p_logger;

    static std::filesystem::path _filePath;
};


} // namespace cie::utils


#endif