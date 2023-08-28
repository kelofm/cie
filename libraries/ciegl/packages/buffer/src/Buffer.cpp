// --- Internal Includes ---
#include "packages/buffer/inc/Buffer.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::gl {


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer,TData>
Buffer<TData,TDataContainer>::Buffer() :
    utils::IDObject<Size>(std::numeric_limits<Size>().max()),
    _drawMode(GL_STATIC_DRAW)
{
    CIE_BEGIN_EXCEPTION_TRACING

    GLuint tempID;
    glGenBuffers(1, &tempID);
    this->setID(tempID);

    CIE_END_EXCEPTION_TRACING
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
Buffer<TData,TDataContainer>::~Buffer()
{
    const GLuint id = this->getID();
    glDeleteBuffers(1, &id);
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
inline void
Buffer<TData,TDataContainer>::resize(Size byteCount)
{
    CIE_BEGIN_EXCEPTION_TRACING

    glNamedBufferData(GLuint(this->_id),
                      byteCount,
                      NULL,
                      _drawMode);

    CIE_END_EXCEPTION_TRACING
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
inline void
Buffer<TData,TDataContainer>::write(ConstDataIterator begin, Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!count) [[unlikely]]
        return;

    glNamedBufferData(
        GLuint(this->_id),
        count * sizeof(data_type),
        &*begin,
        _drawMode
   );

    CIE_END_EXCEPTION_TRACING
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
inline void
Buffer<TData,TDataContainer>::write(Size bufferBegin, ConstDataIterator begin, Size count)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!count) [[unlikely]]
        return;

    glNamedBufferSubData(
        GLuint(this->_id),
        bufferBegin,
        count * sizeof(data_type),
        &*begin
   );

    CIE_END_EXCEPTION_TRACING
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
void
Buffer<TData,TDataContainer>::setDrawMode(Size drawMode)
{
    _drawMode = drawMode;
}


template <class TData, class TDataContainer>
requires concepts::Container<TDataContainer, TData>
const GLuint
Buffer<TData,TDataContainer>::getDrawMode() const
{
    return _drawMode;
}


template class Buffer<float, DynamicArray<float>>;

template class Buffer<float, AttributeContainer>;

template class Buffer<GLint, DynamicArray<GLint>>;

template class Buffer<GLuint, DynamicArray<GLuint>>;


} // namespace cie::gl