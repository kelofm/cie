// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/BufferManager.hpp"


namespace cie::gl {


template <class TBuffer>
struct BufferUtil {};


template <>
struct BufferUtil<VertexBuffer>
{
    using Buffer = VertexBuffer;

    using Iterator = Buffer::ConstDataIterator;

    static void write(Iterator begin, Size count, BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundVertexBuffer)
        r_manager._p_boundVertexBuffer->write(begin, count);
    }

    static void write(Size bufferBegin, Iterator begin, Size count, BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundVertexBuffer)
        r_manager._p_boundVertexBuffer->write(bufferBegin, begin, count);
    }

    static inline Buffer::SharedPointer make(BufferManager& r_manager)
    {
        Buffer::SharedPointer p_buffer = Buffer::SharedPointer(new Buffer);
        p_buffer->setDrawMode(r_manager._drawMode);
        r_manager.registerBuffer(p_buffer, r_manager._vertexBuffers);

        r_manager.logID("create vertex buffer", p_buffer->getID());

        return p_buffer;
    }

    static inline void add(Buffer::SharedPointer p_buffer, BufferManager& r_manager)
    {
        r_manager.registerBuffer(p_buffer, r_manager._vertexBuffers);
    }

    static inline void bind(Buffer::SharedPointer p_buffer, bool log, BufferManager& r_manager)
    {
        if (log) r_manager.logID("Bind vertex buffer", p_buffer->getID());

        const bool exists = std::find(
            r_manager._vertexBuffers.begin(),
            r_manager._vertexBuffers.end(),
            p_buffer
        ) != r_manager._vertexBuffers.end();

        CIE_CHECK(exists, "Attempt to bind buffer that doesn't belong to this manager!")

        glBindBuffer(GL_ARRAY_BUFFER, p_buffer->getID());
        r_manager._p_boundVertexBuffer = p_buffer;
    }

    static bool hasBoundBuffer(const BufferManager& r_manager)
    {return r_manager._p_boundVertexBuffer != nullptr;}

    static inline Buffer::WeakConstPointer getBoundBuffer(const BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundVertexBuffer)
        return r_manager._p_boundVertexBuffer;
    }

    static inline Buffer::WeakPointer getBoundBuffer(BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundVertexBuffer)
        return r_manager._p_boundVertexBuffer;
    }
}; // struct BufferUtil<VertexBuffer>


template <>
struct BufferUtil<ElementBuffer>
{
    using Buffer = ElementBuffer;

    using Iterator = Buffer::ConstDataIterator;

    static void write(Iterator begin, Size count, BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundElementBuffer)
        r_manager._p_boundElementBuffer->write(begin, count);
    }

    static void write(Size bufferBegin, Iterator begin, Size count, BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundElementBuffer)
        r_manager._p_boundElementBuffer->write(bufferBegin, begin, count);
    }

    static inline Buffer::SharedPointer make(BufferManager& r_manager)
    {
        Buffer::SharedPointer p_buffer = Buffer::SharedPointer(new Buffer);
        p_buffer->setDrawMode(r_manager._drawMode);
        r_manager.registerBuffer(p_buffer, r_manager._elementBuffers);

        r_manager.logID("create element buffer", p_buffer->getID());

        return p_buffer;
    }

    static inline void add(Buffer::SharedPointer p_buffer, BufferManager& r_manager)
    {
        r_manager.registerBuffer(p_buffer, r_manager._elementBuffers);
    }

    static inline void bind(Buffer::SharedPointer p_buffer, bool log, BufferManager& r_manager)
    {
        if (log) r_manager.logID("Bind element buffer", p_buffer->getID());

        const bool exists = std::find(
            r_manager._elementBuffers.begin(),
            r_manager._elementBuffers.end(),
            p_buffer
        ) != r_manager._elementBuffers.end();

        CIE_CHECK(exists, "Attempt to bind buffer that doesn't belong to this manager!")

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, p_buffer->getID());
        r_manager._p_boundElementBuffer = p_buffer;
    }

    static bool hasBoundBuffer(const BufferManager& r_manager)
    {return r_manager._p_boundElementBuffer != nullptr;}

    static inline Buffer::WeakConstPointer getBoundBuffer(const BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundElementBuffer)
        return r_manager._p_boundElementBuffer;
    }

    static inline Buffer::WeakPointer getBoundBuffer(BufferManager& r_manager)
    {
        CIE_CHECK_POINTER(r_manager._p_boundElementBuffer)
        return r_manager._p_boundElementBuffer;
    }
}; // struct BufferUtil<ElementBuffer>


BufferManager::BufferManager(utils::Logger& r_logger,
                             const std::string& r_name) :
    utils::Loggee(r_logger, r_name),
    _drawMode(GL_DYNAMIC_DRAW)
{
}


template <class TBuffer>
typename TBuffer::SharedPointer BufferManager::makeBuffer()
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<TBuffer>::make(*this);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Pointer TBufferPtr>
void BufferManager::addBuffer(TBufferPtr p_buffer)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<typename TBufferPtr::element_type>::add(p_buffer, *this);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Pointer TBufferPtr>
void BufferManager::bindBuffer(TBufferPtr p_buffer, bool log)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<typename TBufferPtr::element_type>::bind(p_buffer, log, *this);
    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer>
bool BufferManager::hasBoundBuffer() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<TBuffer>::hasBoundBuffer(*this);
    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer>
typename TBuffer::WeakConstPointer BufferManager::getBoundBuffer() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<TBuffer>::getBoundBuffer(*this);
    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer>
typename TBuffer::WeakPointer BufferManager::getBoundBuffer()
{
    CIE_BEGIN_EXCEPTION_TRACING
    return BufferUtil<TBuffer>::getBoundBuffer(*this);
    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer, concepts::Any<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
void BufferManager::writeToBoundBuffer(TIterator begin,
                                       Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING
    BufferUtil<TBuffer>::write(
        typename TBuffer::ConstDataIterator(begin),
        count,
        *this
    );
    CIE_END_EXCEPTION_TRACING
}


template <class TBuffer, concepts::Any<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
void BufferManager::writeToBoundBuffer(Size bufferBegin,
                                       TIterator begin,
                                       Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING
    BufferUtil<TBuffer>::write(
        bufferBegin,
        typename TBuffer::ConstDataIterator(begin),
        count,
        *this
    );
    CIE_END_EXCEPTION_TRACING
}


template <class TBufferPtr, class TBufferContainer>
void
BufferManager::registerBuffer(TBufferPtr p_buffer, TBufferContainer& r_bufferContainer)
{
    r_bufferContainer.push_back(p_buffer);
}


const typename BufferManager::vertex_buffer_container& BufferManager::vertexBuffers() const
{
    return _vertexBuffers;
}


const typename BufferManager::element_buffer_container& BufferManager::elementBuffers() const
{
    return _elementBuffers;
}


GLuint BufferManager::getDrawMode() const
{
    return _drawMode;
}


void BufferManager::setDrawMode(GLuint drawMode)
{
    _drawMode = drawMode;
    switch (drawMode)
    {
        case GL_STATIC_DRAW:
            log( "Set GL_STATIC_DRAW");
            break;
        case GL_DYNAMIC_DRAW:
            log( "Set GL_DYNAMIC_DRAW");
            break;
        case GL_STREAM_DRAW:
            log( "Set GL_STREAM_DRAW");
            break;
        default:
            log( "Attempt to set invalid draw mode: " + std::to_string(drawMode), LOG_TYPE_ERROR );
            break;
    }
}


// Explicit instantiations for VertexBuffer
template VertexBuffer::SharedPointer BufferManager::makeBuffer<VertexBuffer>();

template void BufferManager::addBuffer<VertexBuffer::SharedPointer>(VertexBuffer::SharedPointer);

template void BufferManager::bindBuffer<VertexBuffer::SharedPointer>(VertexBuffer::SharedPointer,bool);

template bool BufferManager::hasBoundBuffer<VertexBuffer>() const;

template VertexBuffer::WeakPointer BufferManager::getBoundBuffer<VertexBuffer>();

template VertexBuffer::WeakConstPointer BufferManager::getBoundBuffer<VertexBuffer>() const;

template void BufferManager::writeToBoundBuffer<VertexBuffer,typename VertexBuffer::ConstDataIterator>(VertexBuffer::ConstDataIterator,Size);

template void BufferManager::writeToBoundBuffer<VertexBuffer,typename VertexBuffer::ConstDataIterator>(Size,VertexBuffer::ConstDataIterator,Size);

template void BufferManager::writeToBoundBuffer<VertexBuffer,typename VertexBuffer::DataIterator>(VertexBuffer::DataIterator,Size);

template void BufferManager::writeToBoundBuffer<VertexBuffer,typename VertexBuffer::DataIterator>(Size,VertexBuffer::DataIterator,Size);

template void BufferManager::registerBuffer(VertexBuffer::SharedPointer,BufferManager::vertex_buffer_container&);


// Explicit instantiations for ElementBuffer
template ElementBuffer::SharedPointer BufferManager::makeBuffer<ElementBuffer>();

template void BufferManager::addBuffer<ElementBuffer::SharedPointer>(ElementBuffer::SharedPointer);

template void BufferManager::bindBuffer<ElementBuffer::SharedPointer>(ElementBuffer::SharedPointer,bool);

template bool BufferManager::hasBoundBuffer<ElementBuffer>() const;

template ElementBuffer::WeakPointer BufferManager::getBoundBuffer<ElementBuffer>();

template ElementBuffer::WeakConstPointer BufferManager::getBoundBuffer<ElementBuffer>() const;

template void BufferManager::writeToBoundBuffer<ElementBuffer,typename ElementBuffer::ConstDataIterator>(ElementBuffer::ConstDataIterator,Size);

template void BufferManager::writeToBoundBuffer<ElementBuffer,typename ElementBuffer::ConstDataIterator>(Size,ElementBuffer::ConstDataIterator,Size);

template void BufferManager::writeToBoundBuffer<ElementBuffer,typename ElementBuffer::DataIterator>(ElementBuffer::DataIterator,Size);

template void BufferManager::writeToBoundBuffer<ElementBuffer,typename ElementBuffer::DataIterator>(Size,ElementBuffer::DataIterator,Size);

template void BufferManager::registerBuffer(ElementBuffer::SharedPointer,BufferManager::element_buffer_container&);


} // namespace cie::gl