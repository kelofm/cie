// --- Utility Includes ---
#include "packages/ranges/inc/TransformView.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Incldues ---
#include "packages/scene/inc/PartScene.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::gl {


PartScene::PartScene( const std::string& r_name,
                      PartScene::part_container&& r_parts,
                      Shader::SharedPointer p_vertexShader,
                      Shader::SharedPointer p_geometryShader,
                      Shader::SharedPointer p_fragmentShader,
                      utils::Logger& r_logger,
                      VertexBuffer::SharedPointer p_vertexBuffer,
                      ElementBuffer::SharedPointer p_elementBuffer ) :
    Scene( r_name,
           p_vertexShader,
           p_geometryShader,
           p_fragmentShader,
           r_logger,
           p_vertexBuffer,
           p_elementBuffer ),
    _parts( r_parts ),
    _partsRequireUpdate( true )
{
    CIE_BEGIN_EXCEPTION_TRACING

    glEnable( GL_DEPTH_TEST );

    CIE_END_EXCEPTION_TRACING
}


void PartScene::addPart( PartPtr p_part )
{
    this->_parts.push_back( p_part );
    this->_partsRequireUpdate = true;
}


void PartScene::removePart( PartPtr p_part )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto itp_part = std::find( this->_parts.begin(),
                               this->_parts.end(),
                               p_part );

    if ( itp_part != this->_parts.end() )
    {
        this->_parts.erase( itp_part );
        this->_partsRequireUpdate = true;
    }
    else
        CIE_THROW( Exception, "Could not find specified part in container" )

    CIE_END_EXCEPTION_TRACING
}


void PartScene::updateParts()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Allocate memory on the GPU
    Size attributeByteCount = 0;
    Size indexByteCount     = 0;
    Size triangleCount      = 0;
    for ( const auto& rp_part : this->_parts )
    {
        attributeByteCount += rp_part->attributeByteCount();
        indexByteCount     += rp_part->indexByteCount();
        triangleCount      += rp_part->numberOfPrimitives();
    }

    if (auto p_buffer = this->_p_bufferManager->template getBoundBuffer<VertexBuffer>().lock())
        p_buffer->resize(attributeByteCount);
    else
        CIE_THROW(Exception, "Attempt to write to a destroyed buffer")

    if (auto p_buffer = this->_p_bufferManager->template getBoundBuffer<ElementBuffer>().lock())
        p_buffer->resize(indexByteCount);
    else
        CIE_THROW(Exception, "Attempt to write to a destroyed buffer")

    // Write to buffers
    Size attributeByteOffset = 0;
    Size indexByteOffset     = 0;
    Size indexOffset         = 0;

    for ( const auto& rp_part : this->_parts )
    {
        this->_p_bufferManager->writeToBoundBuffer<VertexBuffer>(
            attributeByteOffset,
            rp_part->attributes()->begin(),
            rp_part->attributes()->size()
        );

        // Offset vertex indices (without modifying the original ones)
        auto offsetIndexRange = utils::makeTransformView<Part::index_type>(
            rp_part->indices(),
            [indexOffset](auto index) { return index + indexOffset; }
        );

        this->_p_bufferManager->writeToBoundBuffer<ElementBuffer>(
            indexByteOffset,
            offsetIndexRange.begin(),
            offsetIndexRange.size()
        );

        attributeByteOffset += rp_part->attributeByteCount();
        indexByteOffset     += rp_part->indexByteCount();
        indexOffset         += rp_part->numberOfVertices();
    }

    this->_partsRequireUpdate = false;

    CIE_END_EXCEPTION_TRACING
}





} // namespace cie::gl