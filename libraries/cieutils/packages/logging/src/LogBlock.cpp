// --- Internal Includes ---
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::utils {


LogBlock::LogBlock( const std::string& r_name,
                    Logger& r_logger ) :
    _name( r_name ),
    _timerID( r_logger.startTimer() ),
    _r_logger( r_logger )
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_logger.separate().log( "> " + _name ).increaseIndent().flush();

    CIE_END_EXCEPTION_TRACING
}


LogBlock::~LogBlock()
{
    _r_logger.decreaseIndent();
    logElapsed( "< " + _name + " |", false );
    _r_logger.separate().flush();
}


LogBlock& operator<<( LogBlock& r_block, const std::string& r_message )
{
    CIE_BEGIN_EXCEPTION_TRACING
    return r_block.log( r_message );
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils