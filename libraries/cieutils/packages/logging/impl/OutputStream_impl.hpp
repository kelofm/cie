#ifndef CIE_UTILS_OUTPUT_STREAM_IMPL_HPP
#define CIE_UTILS_OUTPUT_STREAM_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie {


template <class TMessage>
requires (concepts::StringStreamable<TMessage> || (concepts::Container<TMessage> && concepts::StringStreamable<typename std::decay<TMessage>::type::value_type>))
OutputStream& operator<<(OutputStream& r_stream, const TMessage& r_message)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_stream.streamInsertOperation(r_message);
    return r_stream;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container TContainer>
requires concepts::StringStreamable<typename TContainer::value_type>
void OutputStream::streamInsertOperation(const TContainer& r_container)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::stringstream stream;
    for (const auto& r_message : r_container)
        stream << r_message << ' ';
    this->streamInsertOperation(stream.str());

    CIE_END_EXCEPTION_TRACING
}


template <concepts::StringStreamable TMessage>
void OutputStream::streamInsertOperation(const TMessage& r_message)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::stringstream stream;
    stream << r_message;
    this->streamInsertOperation(stream.str());

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie

#endif