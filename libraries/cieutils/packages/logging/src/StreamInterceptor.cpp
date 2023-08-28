// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/logging/inc/StreamInterceptor.hpp"


namespace cie {


StreamInterceptor::StreamInterceptor(OutputStream::SharedPointer p_source,
                                     OutputStream::SharedPointer p_sink)
    : OutputStream(), /*creates a dummy stream by default*/
      _p_source(p_source),
      _p_sink(p_sink),
      _p_swapBuffer(p_source->rdbuf())
{
    CIE_BEGIN_EXCEPTION_TRACING

    _p_source->rdbuf(this);

    CIE_END_EXCEPTION_TRACING
}


StreamInterceptor::~StreamInterceptor()
{
    _p_source->rdbuf(_p_swapBuffer);
}


int StreamInterceptor::overflow(int c)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (c != EOF)
    {
        if (this->pbase() == this->epptr())
        {
            char copy = static_cast<char>(c);
            _p_sink->write(&copy, 1);
        }
        else
            this->sputc(static_cast<char>(c));
    }

    return 0;

    CIE_END_EXCEPTION_TRACING
}


std::streamsize StreamInterceptor::xsputn(const char* p_message, std::streamsize messageSize)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _p_sink->write(p_message, messageSize);
    return messageSize;

    CIE_END_EXCEPTION_TRACING
}


int StreamInterceptor::sync()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (this->pbase() != this->pptr())
    {
        this->_p_sink->write(this->pbase(), this->pptr() - this->pbase());
        this->setp(this->pbase(), this->epptr());
    }

    return 0;

    CIE_END_EXCEPTION_TRACING
}





} // namespace cie