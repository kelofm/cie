// --- External Includes ---
#include <glm/glm.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/file/inc/Triangulated3DPart.hpp"

// --- STL Includes ---
#include <fstream>


namespace cie::gl {


void Triangulated3DPart::writeSTL( const std::filesystem::path& r_path ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK( this->dimension() == 3, "writing to STL is only implemented in 3D" )

    // Create file
    std::ofstream file;
    file.open( r_path, std::ios::binary );

    if ( !file.is_open() )
        CIE_THROW( Exception, "Could not create file: " + r_path.string() )

    // Header
    char header = ' ';
    for ( Size i=0; i<80; ++i )
        file.write( reinterpret_cast<char*>(&header), sizeof(header) );

    // Number of triangles
    uint32_t numberOfTriangles = this->numberOfPrimitives();
    file.write( reinterpret_cast<char*>( &numberOfTriangles ), sizeof(numberOfTriangles) );

    // Attributes
    uint16_t attributeByteCount = this->attributeByteCount();
    
    for ( Size triangleIndex=0; triangleIndex<numberOfTriangles; ++triangleIndex )
    {
        Size triangleBegin = triangleIndex * this->primitiveVertexSize();
        
        Size p0Begin       = this->_indices[triangleBegin] * this->vertexAttributeSize();
        Size p1Begin       = this->_indices[triangleBegin + 1] * this->vertexAttributeSize();
        Size p2Begin       = this->_indices[triangleBegin + 2] * this->vertexAttributeSize();

        // Compute normal
        glm::vec3 p0( _p_attributes->at(p0Begin), _p_attributes->at(p0Begin + 1), _p_attributes->at(p0Begin + 2) );

        glm::vec3 normal = glm::normalize(glm::cross(
            glm::vec3( _p_attributes->at(p2Begin), _p_attributes->at(p2Begin + 1), _p_attributes->at(p2Begin + 2) ) - p0,
            glm::vec3( _p_attributes->at(p1Begin), _p_attributes->at(p1Begin + 1), _p_attributes->at(p1Begin + 2) ) - p0
        ));

        // Write normal
        for ( Size dim=0; dim<this->dimension(); ++dim )
            file.write( reinterpret_cast<char*>(&normal[dim]), sizeof(Part::attribute_type) );

        // Write points
        for ( Size pointBegin : { p0Begin, p1Begin, p2Begin } )
            file.write( reinterpret_cast<const char*>( &_p_attributes->at(pointBegin) ), this->vertexAttributeSize() * sizeof( Triangulated3DPart::attribute_type ) );

        // Write attribute byte count
        file.write( reinterpret_cast<char*>( &attributeByteCount ), sizeof( attributeByteCount ) );
    }

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl