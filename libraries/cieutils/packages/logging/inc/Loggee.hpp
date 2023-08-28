#ifndef CIE_CIEUTILS_LOGGEE_HPP
#define CIE_CIEUTILS_LOGGEE_HPP

// --- Internal Includes ---
#include "packages/logging/inc/Logger.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/types/inc/NamedObject.hpp"

// --- STL Includes ---
#include <string>


// Logger message flags
namespace cie {
const size_t LOG_TYPE_MESSAGE   = 0;
const size_t LOG_TYPE_WARNING   = 1;
const size_t LOG_TYPE_ERROR     = 2;
}


namespace cie::utils {


class Loggee : public NamedObject
{
public:
    Loggee(Logger& r_logger,
           const std::string& r_instanceName);

    ~Loggee();

    void log(const std::string& message,
             size_t messageType = LOG_TYPE_MESSAGE);

    void logID(const std::string& message,
               size_t id,
               size_t messageType = LOG_TYPE_MESSAGE);

    [[nodiscard]] size_t tic();

    void toc(size_t timerID,
             bool reset = true);

    void toc(const std::string& message,
             size_t timerID,
             bool reset = true);

    LogBlock newBlock( const std::string& r_name ) { return _logger->newBlock(r_name); }

    void separate();

    template <class T>
    friend Loggee& operator<<(Loggee& r_loggee, const T& r_object)
    {r_loggee.logger() << r_object; return r_loggee;}

protected:
    Logger& logger();

private:
    Logger*         _logger;
    size_t          _timerID;
};


} // namespace cie::utils

#endif