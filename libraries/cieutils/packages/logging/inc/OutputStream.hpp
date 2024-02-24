#ifndef CIE_UTILS_OUTPUT_STREAM_HPP
#define CIE_UTILS_OUTPUT_STREAM_HPP

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/compile_time/packages/concepts/inc/streamable.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <ostream>
#include <streambuf>
#include <memory>


namespace cie {


/// @addtogroup cieutils
class OutputStream
{
private:
    template <class T>
    class FlaggedDeleter
    {
    public:
        FlaggedDeleter() : _shouldDelete(true) {}

        FlaggedDeleter(bool shouldDelete) : _shouldDelete(shouldDelete) {}

        FlaggedDeleter(FlaggedDeleter<T>&& r_rhs) = default;

        FlaggedDeleter(const FlaggedDeleter<T>& r_rhs) = default;

        FlaggedDeleter<T>& operator=(FlaggedDeleter<T>&& r_rhs) = default;

        FlaggedDeleter<T>& operator=(const FlaggedDeleter<T>& r_rhs) = default;

        void operator()(T* p) {if (_shouldDelete) delete p;}

    private:
        bool _shouldDelete;
    }; // struct FlaggedDeleter

    using StreamBufPtr = std::unique_ptr<std::streambuf,FlaggedDeleter<std::streambuf>>;

    using OStreamPtr   = std::unique_ptr<std::ostream,FlaggedDeleter<std::ostream>>;

    using StreamOperator = std::ostream& (*)(std::ostream&);

public:
    CIE_DEFINE_CLASS_POINTERS(OutputStream)

public:
    /// Create a no-op dummy stream
    OutputStream();

    /// Create an output stream using the specified buffer
    OutputStream(std::streambuf* p_buffer);

    /// Create an output stream wrapper around a statically or dynamically allocated stream
    OutputStream(std::ostream* p_stream,
                 bool isDynamic = false);

    OutputStream(OutputStream&& r_rhs) = default;

    OutputStream(const OutputStream& r_rhs) = delete;

    OutputStream& operator=(OutputStream&& r_rhs) = default;

    OutputStream& operator=(const OutputStream& r_rhs) = delete;

    virtual ~OutputStream() {}

    virtual void write(const char* p_message, std::streamsize messageSize);

    virtual void flush();

    /// Get stream buffer of the wrapped stream
    virtual std::streambuf* rdbuf();

    /// Set stream buffer of the wrapped stream
    virtual void rdbuf(std::streambuf* p_buffer);

    friend OutputStream& operator<<(OutputStream& r_stream, StreamOperator op);

    template <class TMessage>
    requires (concepts::StringStreamable<TMessage> || (concepts::Container<TMessage> && concepts::StringStreamable<typename std::decay<TMessage>::type::value_type>))
    friend OutputStream& operator<<(OutputStream& r_stream, const TMessage& r_t);

public:
    /// TODO: ensure prioritized initialization in static functions
    static OutputStream cout;

    /// TODO: ensure prioritized initialization in static functions
    static OutputStream cerr;

protected:
    template <concepts::Container TContainer>
    requires concepts::StringStreamable<typename TContainer::value_type>
    void streamInsertOperation(const TContainer& r_container);

    template <concepts::StringStreamable TMessage>
    void streamInsertOperation(const TMessage& r_message);

    virtual void streamInsertOperation(const std::string& r_message);

    /// Forward all input operations to the underlying std::ostream
    template <class T>
    inline void forwardToWrappedStream(T&& r_object)
    {
        CIE_CHECK_POINTER(_p_stream) // required because static OutputStreams like cout and cerr might not be initialized yet in static functions
        *(_p_stream) << std::forward<T>(r_object);
    }

private:
    struct DummyBuffer : public std::streambuf
    {
        virtual int overflow(int dummy) override
        {return dummy;}

        virtual std::streamsize xsputn(const char*, std::streamsize messageSize) override
        {return messageSize;}
    };

    StreamBufPtr _p_buffer;

    OStreamPtr   _p_stream;
}; // class OutputStream


} // namespace cie::utils

#include "packages/logging/impl/OutputStream_impl.hpp"

#endif