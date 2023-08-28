// --- Internal Includes ---
#include "packages/logging/inc/OutputStream.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie {


OutputStream::OutputStream()
    : _p_buffer(new OutputStream::DummyBuffer,
                OutputStream::FlaggedDeleter<std::streambuf>(true))
{
    _p_stream = OutputStream::OStreamPtr(
        new std::ostream(_p_buffer.get()),
        OutputStream::FlaggedDeleter<std::ostream>(true)
    );
}


OutputStream::OutputStream(std::streambuf* p_buffer)
    : _p_buffer(p_buffer,
                OutputStream::FlaggedDeleter<std::streambuf>(false))
{
    _p_stream = OutputStream::OStreamPtr(
        new std::ostream(_p_buffer.get()),
        OutputStream::FlaggedDeleter<std::ostream>(true)
    );
}


void OutputStream::write(const char* p_message, std::streamsize messageSize)
{
    _p_stream->write(p_message, messageSize);
}


void OutputStream::flush()
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::flush(*_p_stream);

    CIE_END_EXCEPTION_TRACING
}


std::streambuf* OutputStream::rdbuf()
{
    return _p_stream->rdbuf();
}


void OutputStream::rdbuf(std::streambuf* p_buffer)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_buffer = OutputStream::StreamBufPtr(p_buffer,
                                                 OutputStream::FlaggedDeleter<std::streambuf>(false));
    _p_stream->rdbuf(p_buffer);

    CIE_END_EXCEPTION_TRACING
}


OutputStream::OutputStream(std::ostream* p_stream, bool isDynamic)
    : _p_buffer(p_stream->rdbuf(),
                OutputStream::FlaggedDeleter<std::streambuf>(false)),
      _p_stream(p_stream,
                OutputStream::FlaggedDeleter<std::ostream>(isDynamic))
{
}


OutputStream& operator<<(OutputStream& r_stream, OutputStream::StreamOperator op)
{
    r_stream.streamInsertOperation(op);
    return r_stream;
}


void OutputStream::streamInsertOperation(const std::string& r_message)
{
    this->forwardToWrappedStream(r_message);
}


OutputStream OutputStream::cout(&std::cout);


OutputStream OutputStream::cerr(&std::cerr);


} // namespace cie