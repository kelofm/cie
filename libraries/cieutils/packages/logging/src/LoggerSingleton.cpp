// --- Utility Includes ---
#include "packages/logging/inc/StreamInterceptor.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "cmake_variables.hpp"

// --- Internal Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::utils {


const std::filesystem::path detail::DEFAULT_LOGGER_FILE_PATH = getOutputPath() / "cie.log";


// Construct the static logger if a non-empty file path is specified
LoggerPtr constructLoggerSingleton(const std::filesystem::path& r_filePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    LoggerPtr p_logger;
    if (!r_filePath.empty())
    {
        p_logger = LoggerPtr(new Logger(
            r_filePath
        ));

        auto p_loggerStream = std::static_pointer_cast<OutputStream>(p_logger);

        // Reroute std::cout
        auto p_cout = OutputStream::SharedPointer(new OutputStream(&std::cout, false));
        auto p_coutCopy = OutputStream::SharedPointer(new OutputStream(std::cout.rdbuf()));

        auto p_coutRedirector = OutputStream::SharedPointer(new StreamInterceptor(
            p_cout,
            p_loggerStream
        ));

        // Reroute std::cerr
        auto p_cerr = OutputStream::SharedPointer(new OutputStream(&std::cerr, false));

        auto p_cerrRedirector = OutputStream::SharedPointer(new StreamInterceptor(
            p_cerr,
            p_loggerStream
        ));

        // Reroute std::clog
        auto p_clog = OutputStream::SharedPointer(new OutputStream(&std::clog, false));

        auto p_clogRedirector = OutputStream::SharedPointer(new StreamInterceptor(
            p_clog,
            p_loggerStream
        ));

        // Add original std::cout as stream
        p_logger->addStream(p_coutRedirector);
        p_logger->addStream(p_cerrRedirector);
        p_logger->addStream(p_clogRedirector);

        // Keep the original cout
        p_logger->addStream(p_coutCopy);
    }

    return p_logger;

    CIE_END_EXCEPTION_TRACING
}


// Empty initialization
LoggerPtr LoggerSingleton::_p_logger = constructLoggerSingleton("");


std::filesystem::path LoggerSingleton::_filePath = "";


void LoggerSingleton::initialize(const std::filesystem::path& r_filePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!_p_logger)
    {
        _p_logger = constructLoggerSingleton(r_filePath);
        _filePath = r_filePath;
    }
    else if (_filePath != r_filePath)
        CIE_THROW(Exception, "Attempt to reinitialize LoggerSingleton")

    CIE_END_EXCEPTION_TRACING
}


Logger& LoggerSingleton::get()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!_p_logger)
        initialize();

    return *_p_logger;

    CIE_END_EXCEPTION_TRACING
}


Logger& LoggerSingleton::get(const std::filesystem::path& r_filePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!_p_logger)
        initialize(r_filePath);

    return *_p_logger;

    CIE_END_EXCEPTION_TRACING
}


LoggerPtr LoggerSingleton::getPtr()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!_p_logger)
        initialize();

    return _p_logger;

    CIE_END_EXCEPTION_TRACING
}


LoggerPtr LoggerSingleton::getPtr(const std::filesystem::path& r_filePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!_p_logger)
        initialize(r_filePath);

    return _p_logger;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils