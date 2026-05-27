#pragma once

// --- Utility Includes ---
#include "packages/logging/inc/Logger.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <string>


namespace cie::utils {


/// @ingroup cieutils
class LogBlock
{
public:
    LogBlock() = delete;

    LogBlock(
        const std::string& rName,
        Logger& rLogger);

    ~LogBlock();

    LogBlock& log(const std::string& rMessage);

    LogBlock& warn(const std::string& rMessage);

    LogBlock& error(const std::string& rMessage);

    LogBlock& logElapsed(
        const std::string& rMessage,
        bool reset = true,
        OptionalRef<const Color> rMaybeColor = {});

protected:
    std::string _name;
    Size        _timerID;
    Logger&     _rLogger;
};


LogBlock& operator<<( LogBlock& r_block, const std::string& rMessage );


} // namespace cie::utils
