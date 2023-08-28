#ifndef CIE_CIEUTILS_EXCEPTIONS_IMPL_HPP
#define CIE_CIEUTILS_EXCEPTIONS_IMPL_HPP

// --- STL Includes ---
#include <sstream>


namespace cie {


template <concepts::CIEException TException>
TException exceptionFactory(const String& r_location, const String& r_message, const String&) noexcept
{return TException(r_location, r_message);}


template <concepts::STLException TException>
TException exceptionFactory(const String& r_location, const String& r_message, const String& r_name) noexcept
{return TException((r_location + "\nUnhandled " + r_name + ": " + r_message));}


template <concepts::CIEException TException>
TException exceptionFactory(TException& r_exception, const String& r_location, const String& r_additionalMessage) noexcept
{return TException(r_exception, r_location, r_additionalMessage);}


template <concepts::STLException TException>
TException exceptionFactory(TException& r_exception, const String&, const String& r_additionalMessage) noexcept
{
    std::stringstream stream;
    stream << r_exception.what() << std::endl << r_additionalMessage;
    return TException(stream.str());
}


} // namespace cie


#endif