#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <stdexcept>
#include <string>
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

    Exception(const Exception& r_base, const String& rLocation, const String& r_additionalMessage);

    Exception(const String& location, const String& message);

    Exception(const String& r_what, Size stackLevel);

    Exception(const String& r_what);

    Exception(Exception&& rRHS) = default;

    Exception(const Exception& rRHS);

    Exception& operator=(Exception&& rRHS) = delete;

    Exception& operator=(const Exception& rRHS) = delete;

    virtual ~Exception() {}

    const char* what() const noexcept override;

    virtual std::string name() const;

    Size stackLevel() const;

    template <concepts::CIEException TException>
    friend TException exceptionFactory(const String& rLocation, const String& rMessage, const String& rName) noexcept;

    template <concepts::CIEException TException>
    friend TException exceptionFactory(TException& rException, const String& rLocation, const String& r_additionalMessage) noexcept;

private:
    void setStackLevel(Size stackLevel);

    void setMessage(std::string&& rMessage);

private:
    Size _stackLevel;

    String _what;
};


/// @ingroup cieutils
struct NullPtrException : public Exception
{
    NullPtrException(const String& rLocation, const String& rMessage);

    NullPtrException(NullPtrException&& rRHS) = default;

    NullPtrException(const NullPtrException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct AbstractCallException : public Exception
{
    AbstractCallException(const String& rLocation, const String& rMessage);

    AbstractCallException(AbstractCallException&& rRHS) = default;

    AbstractCallException(const AbstractCallException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct NotImplementedException : public Exception
{
    NotImplementedException(const String& rLocation, const String& rMessage);

    NotImplementedException(NotImplementedException&& rRHS) = default;

    NotImplementedException(const NotImplementedException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct OutOfRangeException : public Exception
{
    OutOfRangeException(const String& rLocation, const String& rMessage);

    OutOfRangeException(OutOfRangeException&& rRHS) = default;

    OutOfRangeException(const OutOfRangeException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct DivisionByZeroException : public Exception
{
    DivisionByZeroException(const String& rLocation, const String& rMessage);

    DivisionByZeroException(DivisionByZeroException&& rRHS) = default;

    DivisionByZeroException(const DivisionByZeroException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct GeometryException : public Exception
{
    GeometryException(const String& rLocation, const String& rMessage);

    GeometryException(GeometryException&& rRHS) = default;

    GeometryException(const GeometryException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
struct MemoryException : public Exception
{
    MemoryException(const String& rLocation, const String& rMessage);

    MemoryException(MemoryException&& rRHS) = default;

    MemoryException(const MemoryException& rRHS) = default;

    virtual std::string name() const override;
};


/// @ingroup cieutils
template <concepts::STLException TException>
TException exceptionFactory(
    const String& rLocation,
    const String& rMessage,
    const String& rName) noexcept;


/// @ingroup cieutils
template <concepts::STLException TException>
TException exceptionFactory(
    TException& rException,
    const String& rLocation,
    const String& rAdditionalMessage) noexcept;


} // namespace cie

#include "packages/exceptions/impl/exceptions_impl.hpp"
