#ifndef CIE_UTILS_STREAM_INTERCEPTOR_HPP
#define CIE_UTILS_STREAM_INTERCEPTOR_HPP

// --- Internal Includes ---
#include "packages/logging/inc/OutputStream.hpp"


namespace cie {


class StreamInterceptor : public std::streambuf, public OutputStream
{
public:
    StreamInterceptor() = delete;

    StreamInterceptor(OutputStream::SharedPointer p_source,
                      OutputStream::SharedPointer p_sink);

    StreamInterceptor(StreamInterceptor&& r_rhs) = default;

    StreamInterceptor(const StreamInterceptor& r_rhs) = delete;

    StreamInterceptor& operator=(StreamInterceptor&& r_rhs) = delete;

    StreamInterceptor& operator=(const StreamInterceptor& r_rhs) = delete;

    virtual ~StreamInterceptor() override;

protected:
    virtual int overflow(int c) override;

    virtual std::streamsize xsputn(const char* p_message, std::streamsize messageSize) override;

    virtual int sync() override;

private:
    OutputStream::SharedPointer _p_source;
    OutputStream::SharedPointer _p_sink;
    std::streambuf*             _p_swapBuffer;
}; // class StreamInterceptor


} // namespace cie


#endif