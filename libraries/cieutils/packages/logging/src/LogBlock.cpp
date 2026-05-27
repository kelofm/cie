// --- Utility Includes ---
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/types/inc/Color.hpp"

// --- STL Includes ---
#include <format>


namespace cie::utils {


LogBlock::LogBlock(
    const std::string& rName,
    Logger& rLogger)
        :   _name(rName),
            _timerID(rLogger.startTimer()),
            _rLogger(rLogger) {
                CIE_BEGIN_EXCEPTION_TRACING
                    rLogger
                        .log(std::format("> {}", _name), RGBAColor::TUMGray)
                        .increaseIndent()
                        .flush();
                CIE_END_EXCEPTION_TRACING
}


LogBlock::~LogBlock() {
    _rLogger.decreaseIndent();
    logElapsed(
        std::format("< {} |", _name),
        false,
        RGBAColor::TUMGray);
    _rLogger.flush();
}


LogBlock& operator<<(
    LogBlock& rBlock,
    const std::string& rMessage) {
        CIE_BEGIN_EXCEPTION_TRACING
            return rBlock.log(rMessage);
        CIE_END_EXCEPTION_TRACING
}


LogBlock& LogBlock::log(const std::string& rMessage)
{
    _rLogger.log(rMessage);
    return *this;
}


LogBlock& LogBlock::warn(const std::string& rMessage)
{
    _rLogger.warn(rMessage);
    return *this;
}


LogBlock& LogBlock::error(const std::string& rMessage)
{
    _rLogger.error(rMessage);
    return *this;
}


LogBlock& LogBlock::logElapsed(
    const std::string& rMessage,
    bool reset,
    OptionalRef<const Color> rMaybeColor) {
        _rLogger.logElapsed(
            rMessage,
            _timerID,
            reset,
            rMaybeColor);
        return *this;
}


} // namespace cie::utils
