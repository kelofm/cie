#ifndef CIE_CIEUTILS_EXCEPTIONS_HPP
#define CIE_CIEUTILS_EXCEPTIONS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <stdexcept>
#include <string>
#include <deque>
#include <memory>


namespace cie {
class Exception;
} // namespace cie


namespace cie::concepts {

/// @ingroup cieutils
template <class T>
concept CIEException
= std::derived_from<T,cie::Exception>;

/// @ingroup cieutils
template <class T>
concept STLException
= std::derived_from<T,std::exception>
  && !CIEException<T>;

} // namespace cie::concepts


namespace cie {


/// @ingroup cieutils
class Exception : public std::exception
{
public:
    Exception() = delete;

    Exception(const String& location, const String& message, Size stackLevel);

    Exception(const Exception& r_base, const String& r_location, const String& r_additionalMessage);

    Exception(const String& location, const String& message);

    Exception(const String& r_what, Size stackLevel);

    Exception(const String& r_what);

    Exception(Exception&& r_rhs) = default;

    Exception(const Exception& r_rhs);

    Exception& operator=(Exception&& r_rhs) = delete;

    Exception& operator=(const Exception& r_rhs) = delete;

    virtual ~Exception() {}

    const char* what() const noexcept override;

    virtual std::string name() const;

    Size stackLevel() const;

    template <concepts::CIEException TException>
    friend TException exceptionFactory(const String& r_location, const String& r_message, const String& r_name) noexcept;

    template <concepts::CIEException TException>
    friend TException exceptionFactory(TException& r_exception, const String& r_location, const String& r_additionalMessage) noexcept;

private:
    void setStackLevel(Size stackLevel);

    void setMessage(std::string&& r_message);

private:
    Size _stackLevel;

    String _what;
};


/// @ingroup cieutils
struct NullPtrException : public Exception
{
    NullPtrException(const String& r_location, const String& r_message);

    NullPtrException(NullPtrException&& r_rhs) = default;

    NullPtrException(const NullPtrException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct AbstractCallException : public Exception
{
    AbstractCallException(const String& r_location, const String& r_message);

    AbstractCallException(AbstractCallException&& r_rhs) = default;

    AbstractCallException(const AbstractCallException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct NotImplementedException : public Exception
{
    NotImplementedException(const String& r_location, const String& r_message);

    NotImplementedException(NotImplementedException&& r_rhs) = default;

    NotImplementedException(const NotImplementedException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct OutOfRangeException : public Exception
{
    OutOfRangeException(const String& r_location, const String& r_message);

    OutOfRangeException(OutOfRangeException&& r_rhs) = default;

    OutOfRangeException(const OutOfRangeException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct DivisionByZeroException : public Exception
{
    DivisionByZeroException(const String& r_location, const String& r_message);

    DivisionByZeroException(DivisionByZeroException&& r_rhs) = default;

    DivisionByZeroException(const DivisionByZeroException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct GeometryException : public Exception
{
    GeometryException(const String& r_location, const String& r_message);

    GeometryException(GeometryException&& r_rhs) = default;

    GeometryException(const GeometryException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct MemoryException : public Exception
{
    MemoryException(const String& r_location, const String& r_message);

    MemoryException(MemoryException&& r_rhs) = default;

    MemoryException(const MemoryException& r_rhs) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
template <concepts::STLException TException>
TException exceptionFactory(const String& r_location, const String& r_message, const String& r_name) noexcept;


/// @ingroup cieutils
template <concepts::STLException TException>
TException exceptionFactory(TException& r_exception, const String& r_location, const String& r_additionalMessage) noexcept;


} // namespace cie

#include "packages/exceptions/impl/exceptions_impl.hpp"

#endif