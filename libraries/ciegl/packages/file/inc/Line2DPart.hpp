#ifndef CIE_GL_FILE_LINE_2D_PART_HPP
#define CIE_GL_FILE_LINE_2D_PART_HPP

// --- Internal Includes ---
#include "packages/file/inc/Part.hpp"


namespace cie::gl {


class Line2DPart : public Part
{
public:
    virtual void load( std::istream& r_stream ) override;

    virtual Size dimension() const override;

    virtual Size vertexAttributeSize() const override;

    virtual Size primitiveVertexSize() const override;

    virtual void writeSTL( const std::filesystem::path& r_path ) const override;
};


} // namespace cie::gl


#endif