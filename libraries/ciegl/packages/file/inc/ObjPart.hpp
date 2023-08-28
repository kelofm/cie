#ifndef CIE_GL_OBJ_PART_HPP
#define CIE_GL_OBJ_PART_HPP

// --- Internal Includes ---
#include "packages/file/inc/Triangulated3DPart.hpp"


namespace cie::gl {


/**
 * Load geometry from an .obj file
 * 
 * Note: vertex positions are the only attributes loaded
 */
class ObjPart final : public Triangulated3DPart
{
public:
    ObjPart( std::istream& r_stream );
    ObjPart( const std::filesystem::path& r_filePath );

    void load( std::istream& r_stream ) override;
};


} // namespace cie::gl


#endif