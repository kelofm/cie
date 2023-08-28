#ifndef CIE_UTILS_LOGGING_LOG_BLOCK_IMPL_HPP
#define CIE_UTILS_LOGGING_LOG_BLOCK_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils {


inline LogBlock& LogBlock::log( const std::string& r_message )
{
    _r_logger.log( r_message );
    return *this;
}


inline LogBlock& LogBlock::warn( const std::string& r_message )
{
    _r_logger.warn( r_message );
    return *this;
}


inline LogBlock& LogBlock::error( const std::string& r_message )
{
    _r_logger.error( r_message );
    return *this;
}


inline LogBlock& LogBlock::logElapsed(  const std::string& r_message,
                                        bool reset )
{
    _r_logger.logElapsed(   r_message,
                            _timerID,
                            reset );
    return *this;
}


} // namespace cie::utils


#endif