// --- External Includes ---
#include <glm/glm.hpp>

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/file/inc/STLPart.hpp"

// --- STL Includes ---
#include <cstdint>
#include <vector>
#include <algorithm>
#include <fstream>
#include <numeric>


namespace cie::gl {


STLPart::STLPart( std::istream& r_stream )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->load( r_stream );

    CIE_END_EXCEPTION_TRACING
}


STLPart::STLPart( const std::filesystem::path& r_filePath )
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::ifstream file;
    file.open( r_filePath,
               std::ios::binary );

    if ( !file )
        CIE_THROW( Exception, "Failed to open file: " + r_filePath.string() )

    this->load( file );

    CIE_END_EXCEPTION_TRACING
}


void STLPart::load( std::istream& r_stream )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_attributes->clear();
    this->_indices.clear();

    // Header - 80 chars
    {
        char tmp;

        for ( Size i=0; i<80; ++i )
        {
            r_stream.read( reinterpret_cast<char*>(&tmp), sizeof(tmp) );
        }
    }

    // Number of triangles
    Size numberOfTriangles;

    {
        uint32_t tmp;
        r_stream.read( reinterpret_cast<char*>(&tmp), sizeof(tmp) );

        numberOfTriangles = Size(tmp);
    }

    // Allocate memory for the data
    this->_p_attributes->reserve( numberOfTriangles * this->primitiveAttributeSize() );

    // Triangles and their normals
    using attribute_type = typename STLPart::attribute_type;
    attribute_type tmp;
    uint16_t attributeByteCount;

    for ( Size triangleIndex=0; triangleIndex<numberOfTriangles; ++triangleIndex )
    {
        // Read normal
        for ( Size componentIndex=0; componentIndex<this->dimension(); ++componentIndex )
            r_stream.read( reinterpret_cast<char*>(&tmp), sizeof(tmp) );

        // Read points
        for ( Size pointIndex=0; pointIndex<3; pointIndex++ )
        {
            for ( Size componentIndex=0; componentIndex<this->dimension(); ++componentIndex )
            {
                r_stream.read( reinterpret_cast<char*>(&tmp), sizeof(tmp) );
                this->_p_attributes->push_back( tmp );
            }
        }

        // Read "attribute byte count"
        r_stream.read( reinterpret_cast<char*>(&attributeByteCount), sizeof(attributeByteCount) );
    }

    // Contiguous triangles
    this->_indices.resize( numberOfTriangles * this->primitiveVertexSize() );
    std::iota( this->_indices.begin(),
               this->_indices.end(),
               0 );

    if ( r_stream.peek() != EOF )
        CIE_THROW( Exception, "Reading STL file failed" )

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl