// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/logging/inc/Logger.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>


namespace cie::utils {


namespace detail {

Time getTime()
{
    CIE_BEGIN_EXCEPTION_TRACING

    return std::chrono::steady_clock::now();

    CIE_END_EXCEPTION_TRACING
}

std::string getDate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    return std::string(std::ctime(&t));

    CIE_END_EXCEPTION_TRACING
}

} // namespace detail


Logger::Logger(const std::filesystem::path& r_filePath) :
    _timeLog( {detail::getTime()} ),
    _prefix( "" ),
    _lineEnd("\n"),
    _forceFlush(true),
    _newLine(true)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Create log file
    auto p_file = std::static_pointer_cast<OutputStream>(
        std::make_shared<FileOutputStream>(r_filePath, std::ios::out)
    );
    addStream(p_file);

    // Log time
    logDate( "Log file created on" );

    CIE_END_EXCEPTION_TRACING
}


Logger::~Logger()
{
    noIndent();
    logElapsed("\nLogger ran for", 0, false);
    logDate("Log closed on");
    flush();
}


void Logger::write(const char* p_message, std::streamsize messageSize)
{
    CIE_BEGIN_EXCEPTION_TRACING

    const char* it_begin            = p_message;
    const char* const it_messageEnd = p_message + messageSize;
    const char* it_end              = it_messageEnd;

    while (it_begin < it_messageEnd)
    {
        it_end = std::find(
            it_begin,
            it_messageEnd,
            '\n'
        );

        this->preWrite();

        this->directWriteToAll(it_begin,
                               std::distance(it_begin, it_end));

        if (it_end != it_messageEnd)
            this->postWrite();

        it_begin = it_end + 1;
    }

    if (_forceFlush)
        this->flush();

    CIE_END_EXCEPTION_TRACING
}


LogBlock Logger::newBlock( const std::string& r_name )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return LogBlock(r_name, *this);

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::addStream(OutputStream::SharedPointer p_stream)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _streams.push_back(p_stream);
    return *this;

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::removeStream(OutputStream::SharedPointer p_stream)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it = std::find(_streams.begin(),
                        _streams.end(),
                        p_stream);

    if (it != _streams.end() && it!=_streams.begin())
        _streams.erase(it);

    return *this;

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::forceFlush( bool use )
{
    _forceFlush = use;
    return *this;
}


Logger& Logger::log(const std::string& message)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->write(message.c_str(),
                message.size());
    this->postWrite();
    return *this;

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::warn( const std::string& message )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return log( "WARNING: " + message );

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::logDate( const std::string& message )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return log( message + " " + detail::getDate() );

    CIE_END_EXCEPTION_TRACING
}


size_t Logger::startTimer()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find an idle slot
    size_t slotID = std::distance(_timeLog.begin(),
                                  std::find(_timeLog.begin()+1,
                                            _timeLog.end(),
                                            _timeLog[0]));
    if (slotID == _timeLog.size())
    {
        _timeLog.push_back( _timeLog[0] );
        slotID = _timeLog.size() - 1;
    }

    // Start timer
    _timeLog[slotID] = detail::getTime();
    return slotID;

    CIE_END_EXCEPTION_TRACING
}


size_t Logger::elapsed( size_t slotID, bool reset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check if valid slot
    if ( slotID >= _timeLog.size() )
        error( "Invalid timerID " + std::to_string(slotID) );

    // Compute elapsed time
    size_t t = (size_t)std::chrono::duration_cast<std::chrono::microseconds>(detail::getTime() - _timeLog[slotID]).count();

    // Reset timer if requested
    if (reset)
        _timeLog[slotID] = _timeLog[0];

    return t;

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::logElapsed( const std::string& message,
                            size_t timeID,
                            bool reset )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto dt = elapsed( timeID, reset );
    std::string unit = " [us] ";

    if ( dt > 100000 )
    {
        dt      /= 1000;
        unit    = " [ms] ";

        if ( dt > 100000 )
        {
            dt      /= 1000;
            unit    = " [s] ";

            if ( dt > 6000 )
            {
                dt      /= 60;
                unit    = " [min] ";
            }
        }
    }

    return log( message + " " + std::to_string( dt ) + unit );

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::separate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // TODO: separator behaviour is hard-coded
    const char separatorCharacter   = '-';
    const Size lineWidth            = 60;

    Size numberOfSeparators         = 0;
    if ( lineWidth > _prefix.size() )
        numberOfSeparators = lineWidth - _prefix.size();

    return log(std::string(numberOfSeparators, separatorCharacter));

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::increaseIndent()
{
    _prefix += "|   ";
    return *this;
}


Logger& Logger::decreaseIndent()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( _prefix.size() >= 4 )
        _prefix.resize( _prefix.size() - 4 );
    else
        _prefix.clear();
    return *this;

    CIE_END_EXCEPTION_TRACING
}


Logger& Logger::noIndent()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( !_prefix.empty())
        while( _prefix.back() == '\t' )
            _prefix.pop_back();

    return *this;

    CIE_END_EXCEPTION_TRACING
}


void Logger::streamInsertOperation(const std::string& r_message)
{
    this->log(r_message);
}


void Logger::flush()
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (auto& rp_stream : _streams)
        rp_stream->flush();

    CIE_END_EXCEPTION_TRACING
}


void Logger::preWrite()
{
    if (_newLine)
    {
        this->directWriteToAll(_prefix.c_str(),
                               _prefix.size());
        _newLine = false;
    }
}

void Logger::postWrite()
{
    this->directWriteToAll(_lineEnd.c_str(),
                           _lineEnd.size());
    _newLine = true;
}


Logger& Logger::directWriteToAll(const char* p_message,
                                 std::streamsize messageSize)
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (auto& rp_stream : _streams)
        rp_stream->write(p_message, messageSize);

    return *this;

    CIE_END_EXCEPTION_TRACING
}


Logger& operator<<( Logger& r_logger, const std::string& r_message )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return r_logger.log( r_message );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils