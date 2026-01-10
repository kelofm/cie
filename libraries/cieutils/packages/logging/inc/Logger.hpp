#pragma once

// --- Internal Includes ---
#include "packages/logging/inc/FileOutputStream.hpp"

// --- STL Includes ---
#include <fstream>
#include <iosfwd>
#include <deque>
#include <memory>
#include <string>
#include <chrono>
#include <ctime>
#include <filesystem>


namespace cie::utils {


// Clock definition
namespace detail {
using Time = std::chrono::steady_clock::time_point;

Time getTime();

std::string getDate();
} // namespace detail


class LogBlock;


/// @ingroup cieutils
class Logger : public OutputStream
{
public:
    friend class LogBlock;

public:
    Logger() = delete;
    Logger(const std::filesystem::path& r_filePath);
    Logger(const Logger& copy) = delete;
    virtual ~Logger();

    virtual void write(const char* p_message, std::streamsize messageSize) override;

    [[nodiscard]] LogBlock newBlock( const std::string& r_name );

    Logger& addStream(OutputStream::SharedPointer p_stream);
    Logger& removeStream(OutputStream::SharedPointer p_stream);

    Logger& forceFlush(bool use);

    Logger& log(const std::string& message);

    template <class ...Args>
    Logger& logs(Args&&... args);

    Logger& warn(const std::string& message);

    template <class ExceptionType = std::runtime_error>
    Logger& error(const std::string& r_message);

    Logger& logDate(const std::string& message);
    [[nodiscard]] size_t startTimer();
    size_t elapsed(size_t timeID,
                   bool reset=true);
    Logger& logElapsed(const std::string& message,
                       size_t timeID,
                       bool reset=true);

    Logger& separate();

    Logger& increaseIndent();
    Logger& decreaseIndent();
    Logger& noIndent();

protected:
    virtual void streamInsertOperation(const std::string& r_message) override;

    virtual void flush() override;

    void preWrite();

    void postWrite();

    Logger& directWriteToAll(const char* p_message,
                             std::streamsize messageSize);

    std::deque<OutputStream::SharedPointer> _streams;
    std::deque<detail::Time>                _timeLog;

    std::string                             _prefix;
    std::string                             _lineEnd;
    bool                                    _forceFlush;

private:
    bool                                    _newLine;
};


using LoggerPtr = std::shared_ptr<Logger>;


} // namespace cie::utils


#include "packages/logging/impl/Logger_impl.hpp"
