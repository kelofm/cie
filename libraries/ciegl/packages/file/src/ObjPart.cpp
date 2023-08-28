// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/file/inc/ObjPart.hpp"

// --- STL Includes ---
#include <deque>
#include <array>
#include <sstream>
#include <fstream>
#include <algorithm>


namespace cie::gl {


ObjPart::ObjPart( std::istream& r_stream )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->load( r_stream );

    CIE_END_EXCEPTION_TRACING
}


ObjPart::ObjPart( const std::filesystem::path& r_filePath )
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


void ObjPart::load( std::istream& r_stream )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_attributes->clear();
    this->_indices.clear();

    std::string line, token;
    bool isVertex, isFace;

    // Read next line in file and decide whether to continue parsing it
    // -> either vertex pos data, or vertex connectivity data for triangles
    auto readLine = [&]() -> bool
    {
        isVertex = false;
        isFace   = false;

        std::getline( r_stream, line );

        auto delimiterIndex = line.find( ' ' );
        if ( delimiterIndex==0 || delimiterIndex==line.npos )
            return false;

        token = line.substr( 0, delimiterIndex );
        
        if ( token == "v" )
            isVertex = true;
        else if ( token == "f" )
            isFace = true;
        else
            return false;

        line = line.substr( delimiterIndex+1 );

        // Filter texture and normal indices
        delimiterIndex = line.find( '/' );
        while ( delimiterIndex != line.npos )
        {
            Size nextSpace = line.find( ' ', delimiterIndex );
            if ( nextSpace == line.npos )
                nextSpace = line.size();

            std::fill(
                line.begin() + delimiterIndex,
                line.begin() + nextSpace,
                ' '
            );

            delimiterIndex = line.find( '/' );
        }

        return true;
    };


    while ( r_stream.peek() != EOF )
    {
        if ( !readLine() )
            continue;

        if ( isVertex )
        {
            std::stringstream lineStream( line );

            Size vertexBegin = this->_p_attributes->size();
            this->_p_attributes->resize( vertexBegin + this->dimension() );

            for ( Size dim=0; dim<this->dimension(); ++dim )
                lineStream >> this->_p_attributes->at(vertexBegin + dim);
        }
        else if ( isFace )
        {
            std::stringstream lineStream( line );

            std::deque<Size> indices;
            int tmp;

            while( (lineStream >> tmp) )
            {
                if ( tmp >= 0 )
                    indices.push_back( Size(tmp) );
                else
                    CIE_THROW( Exception, "Unsupported indexing in .obj file" )
                    //indices.push_back( vertices.size() + tmp + 1 );
            }

            // Pad indices
            if ( indices.size() % 3 > 0 )
            {
                indices.push_back(indices[0]);
                if ( indices.size() % 3 > 1 )
                    indices.push_back(indices[1]);
            }

            for ( Size offset=0; offset<indices.size() - 2; offset+=2 )
            {
                for ( Size dim=0; dim<this->dimension(); ++dim )
                    this->_indices.push_back( indices[offset+dim] - 1 );
            }
        }
    }

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl