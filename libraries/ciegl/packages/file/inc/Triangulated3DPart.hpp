#ifndef CIE_GL_TRIANGULATED_3D_PART_HPP
#define CIE_GL_TRIANGULATED_3D_PART_HPP

// --- Internal Includes ---
#include "packages/file/inc/Part.hpp"


namespace cie::gl {


class Triangulated3DPart : public Part
{
public:
    Size dimension() const override { return 3; }

    /*! Get number of attributes that make up a primitive
     * @note this override only accounts for vertex position
     */
    virtual Size vertexAttributeSize() const override { return 3; }

    Size primitiveVertexSize() const override { return 3; }

    virtual void writeSTL( const std::filesystem::path& r_path ) const override;
};


} // namespace cie::gl


#endif