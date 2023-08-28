// --- Internal Includes ---
#include "packages/exceptions/inc/exception.hpp"


namespace cie {


Exception::Exception(const String& r_location, const String& r_message, Size stackLevel)
    : _stackLevel(stackLevel)
{
    _what =
        "#" + std::to_string(_stackLevel) + " "
        + r_location + "\n";

    if (stackLevel == 0)
        _what += "Unhandled " + this->name() + ": ";

    if (!r_message.empty())
        _what += r_message;
}


Exception::Exception(const Exception& r_base, const String& r_location, const String& r_additionalMessage)
    : Exception(r_location, String(r_base.what()) + "\n" + r_additionalMessage, r_base.stackLevel() + 1)
{
}


Exception::Exception(const String& r_location, const String& r_message)
    : Exception(r_location, r_message, 0)
{
}


Exception::Exception(const String& r_what, Size stackLevel)
    : _stackLevel(stackLevel),
      _what(r_what)
{
}


Exception::Exception(const String& r_what)
    : Exception(r_what, 0)
{
}


Exception::Exception(const Exception& r_rhs)
    : Exception(r_rhs._what, r_rhs._stackLevel)
{
}


std::string Exception::name() const
{
    return "Exception";
}


const char* Exception::what() const noexcept
{
    return _what.c_str();
}


Size Exception::stackLevel() const
{
    return _stackLevel;
}


void Exception::setStackLevel(Size stackLevel)
{
    _stackLevel = stackLevel;
}


void Exception::setMessage(String&& r_message)
{
    _what = std::move(r_message);
}


NullPtrException::NullPtrException(const String& r_location, const String& r_message)
    : Exception(r_location, "Attempt to dereference a nullptr:\n" + r_message)
{
}


std::string NullPtrException::name() const
{
    return "NullPtrException";
}


AbstractCallException::AbstractCallException(const String& r_location, const String& r_message)
    : Exception(r_location, "Call to member function of abstact class:\n" + r_message)
{
}


std::string AbstractCallException::name() const
{
    return "AbstractCallException";
}


NotImplementedException::NotImplementedException(const String& r_location, const String& r_message)
    : Exception(r_location, "Call to not implemented function/class:\n" + r_message)
{
}


std::string NotImplementedException::name() const
{
    return "NotImplementedException";
}


OutOfRangeException::OutOfRangeException(const String& r_location,
                                            const String& r_message) :
    Exception(r_location, "Index out-of-range or size mismatch:\n" + r_message)
{
}


std::string OutOfRangeException::name() const
{
    return "OutOfRangeException";
}


DivisionByZeroException::DivisionByZeroException(const String& r_location,
                                                    const String& r_message) :
    Exception(r_location, "Division by zero exception:\n" + r_message)
{
}


std::string DivisionByZeroException::name() const
{
    return "DivisionByZeroException";
}


GeometryException::GeometryException(const String& r_location,
                                        const String& r_message) :
    Exception(r_location, "Geometry exception\n" + r_message)
{
}


std::string GeometryException::name() const
{
    return "GeometryException";
}


MemoryException::MemoryException(const String& r_location,
                                  const String& r_message) :
    Exception(r_location, "Memory exception\n" + r_message)
{
}


std::string MemoryException::name() const
{
    return "MemoryException";
}


} // namespace cie