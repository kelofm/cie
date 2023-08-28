#ifndef CIE_GL_ABS_BUFFER_HPP
#define CIE_GL_ABS_BUFFER_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/types/inc/IDObject.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/AttributeContainer.hpp"


namespace cie::gl {


template <class TData, class TDataContainer = DynamicArray<TData>>
requires concepts::Container<TDataContainer,TData>
class Buffer : public utils::IDObject<Size>
{
public:
    using data_type = TData;

    using data_container_type = TDataContainer;   // <-- contiguous data container

    using DataIterator = typename data_container_type::iterator;

    using ConstDataIterator = typename data_container_type::const_iterator;

    CIE_DEFINE_CLASS_POINTERS(Buffer)

    friend class BufferManager;

public:
    Buffer();

    ~Buffer();

//protected:
    /// Allocate the requested number of bytes
    void resize(Size byteCount);

    /// Resize and write data to the GPU buffer.
    void write(ConstDataIterator begin, Size count);

    /// Resize and write data to the GPU buffer.
    void write(Size bufferBegin, ConstDataIterator begin, Size count);

    void setDrawMode(Size drawMode);

    const GLuint getDrawMode() const;

private:
    Buffer(Buffer&& r_rhs) = delete;

    Buffer(const Buffer& r_rhs) = delete;

    Buffer& operator=(Buffer&& r_rhs) = delete;

    Buffer& operator=(const Buffer& r_rhs) = delete;

private:
    GLuint _drawMode;
}; // class Buffer


using VertexBuffer       = Buffer<GLfloat,AttributeContainer>;

using ElementBuffer      = Buffer<GLuint>;


} // namespace cie::gl

#endif