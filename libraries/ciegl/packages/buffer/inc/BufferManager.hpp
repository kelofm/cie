#ifndef CIE_GL_ABS_BUFFER_MANAGER_HPP
#define CIE_GL_ABS_BUFFER_MANAGER_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/logical.hpp"
#include "packages/logging/inc/Loggee.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- GL Includes ---
#include "packages/buffer/inc/Buffer.hpp"

// --- STL Includes ---
#include <list>


// Avoid compiler weirdness
namespace cie::concepts::detail {
template <class T, class TExclude>
concept Except
= !std::is_same_v<typename std::decay<T>::type, typename std::decay<TExclude>::type>;
} // namespace cie::concepts::detail


namespace cie::gl {


class BufferManager : public utils::Loggee
{
public:
    using vertex_buffer_container   = std::list<VertexBuffer::SharedPointer>;

    using element_buffer_container  = std::list<ElementBuffer::SharedPointer>;

    CIE_DEFINE_CLASS_POINTERS(BufferManager)

public:
    BufferManager(utils::Logger& r_logger = utils::LoggerSingleton::get(),
                  const std::string& r_name = "BufferManager");

    /// Create a new buffer and register it in the internal container.
    template <class TBuffer>
    [[nodiscard]] typename TBuffer::SharedPointer makeBuffer();

    /// Register an existing buffer
    template <concepts::Pointer TBufferPtr>
    void addBuffer(TBufferPtr p_buffer);

    /// Bind vertex buffer and set bound flag.
    template <concepts::Pointer TBufferPtr>
    void bindBuffer(TBufferPtr p_buffer, bool log = true);

    template <class TBuffer>
    bool hasBoundBuffer() const;

    template <class TBuffer>
    typename TBuffer::WeakPointer getBoundBuffer();

    template <class TBuffer>
    typename TBuffer::WeakConstPointer getBoundBuffer() const;

    /// Write a contiguous array of values to the bound buffer.
    template <class TBuffer, concepts::Any<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
    void writeToBoundBuffer(TIterator begin, Size count);

    /// Write a contiguous array of values to the bound buffer beginning at the specified position.
    template <class TBuffer, concepts::Any<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
    void writeToBoundBuffer(Size bufferBegin,
                            TIterator begin,
                            Size count);

    /// Copy values to a contiguous container, then write them to the bound buffer.
    template <class TBuffer, concepts::None<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
    void writeToBoundBuffer(TIterator begin, Size count);

    /// Copy values to a contiguous container, then write them to the bound buffer beginning at the specified position.
    template <class TBuffer, concepts::None<typename TBuffer::ConstDataIterator,typename TBuffer::DataIterator> TIterator>
    void writeToBoundBuffer(Size bufferBegin, TIterator begin, Size count);

    GLuint getDrawMode() const;

    void setDrawMode(GLuint drawMode);

    const vertex_buffer_container& vertexBuffers() const;

    const element_buffer_container& elementBuffers() const;

private:
    template <class TBufferPtr, class TBufferContainer>
    void registerBuffer(TBufferPtr p_buffer, TBufferContainer& r_bufferContainer);

private:
    VertexBuffer::SharedPointer _p_boundVertexBuffer;

    ElementBuffer::SharedPointer _p_boundElementBuffer;

    vertex_buffer_container _vertexBuffers;

    element_buffer_container _elementBuffers;

    GLuint _drawMode;

    template <typename>
    friend struct BufferUtil;
};


} // namespace cie::gl

#include "packages/buffer/impl/BufferManager_impl.hpp"

#endif