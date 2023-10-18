#ifndef CIE_UTILS_LOGGING_LOG_BLOCK_HPP
#define CIE_UTILS_LOGGING_LOG_BLOCK_HPP

// --- Utility Includes ---
#include "packages/logging/inc/Logger.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <string>


namespace cie::utils {


/// @addtogroup cieutils
class LogBlock
{
public:
    LogBlock() = delete;
    LogBlock(const std::string& r_name,
             Logger& r_logger);
    ~LogBlock();

    LogBlock& log( const std::string& r_message );
    LogBlock& warn( const std::string& r_message );
    LogBlock& error( const std::string& r_message );

    LogBlock& logElapsed(   const std::string& r_message,
                            bool reset = true );

protected:
    std::string _name;
    Size        _timerID;
    Logger&     _r_logger;
};


LogBlock& operator<<( LogBlock& r_block, const std::string& r_message );


} // namespace cie::utils

#include "packages/logging/impl/LogBlock_impl.hpp"

#endif