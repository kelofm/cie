#ifndef CIE_GL_ABS_BUFFER_MANAGER_IMPL_HPP
#define CIE_GL_ABS_BUFFER_MANAGER_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/resize.hpp"

// --- STL Includes ---
#include <iterator>
#include <algorithm>


namespace cie::gl {


template <class TBuffer, concepts::None<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
void BufferManager::writeToBoundBuffer(TIterator begin, Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if constexpr (std::is_same_v<TBuffer,VertexBuffer>)
    {
        typename TBuffer::data_container_type container;
        utils::reserve(container, count);
        std::copy(begin, begin + count, std::back_inserter(container));
        this->writeToBoundBuffer<TBuffer>(container.cbegin(), count);
    }
    else if constexpr (std::is_same_v<TBuffer,ElementBuffer>)
    {
        typename TBuffer::data_container_type container;
        utils::reserve(container, count);
        std::copy(begin, begin + count, std::back_inserter(container));
        this->writeToBoundBuffer<TBuffer>(container.cbegin(), count);
    }
    //else
    //    static_assert(false, "VertexBuffer or ElementBuffer only!");

    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer, concepts::None<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
void BufferManager::writeToBoundBuffer(Size bufferBegin, TIterator begin, Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if constexpr (std::is_same_v<TBuffer,VertexBuffer>)
    {
        typename TBuffer::data_container_type container;
        utils::reserve(container, count);
        std::copy(begin, begin + count, std::back_inserter(container));
        this->writeToBoundBuffer<TBuffer>(bufferBegin, container.cbegin(), count);
    }
    else if constexpr (std::is_same_v<TBuffer,ElementBuffer>)
    {
        typename TBuffer::data_container_type container;
        utils::reserve(container, count);
        std::copy(begin, begin + count, std::back_inserter(container));
        this->writeToBoundBuffer<TBuffer>(bufferBegin, container.cbegin(), count);
    }
//    else
//        static_assert(false, "VertexBuffer or ElementBuffer only!");

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif