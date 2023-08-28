#ifndef CIE_GL_STL_PART_HPP
#define CIE_GL_STL_PART_HPP

// --- Internal Includes ---
#include "packages/file/inc/Triangulated3DPart.hpp"

// --- STL Includes ---
#include <string>


namespace cie::gl {


class STLPart final : public Triangulated3DPart
{
public:
    STLPart( std::istream& r_stream );
    STLPart( const std::filesystem::path& r_filePath );

    void load( std::istream& r_stream ) override;
};


using STLPartPtr = std::shared_ptr<STLPart>;


} // namespace cie::gl


#endif