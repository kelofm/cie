#ifndef CIE_GL_GENERIC_PART_HPP
#define CIE_GL_GENERIC_PART_HPP

// --- Internal Includes ---
#include "packages/file/inc/Part.hpp"


namespace cie::gl {


class GenericPart : public Part
{
public:
    GenericPart(const Size dimension,
                const Size vertexAttributeSize,
                const Size primitiveVertexSize);

    GenericPart::attribute_container::SharedPointer attributes();

    GenericPart::index_container& indices();

    Size dimension() const override;

    Size vertexAttributeSize() const override;

    Size primitiveVertexSize() const override;

    void writeSTL( const std::filesystem::path& r_path ) const override;

protected:
    const Size _dimension;

    const Size _vertexAttributeSize;

    const Size _primitiveVertexSize;
};


using GenericPartPtr = std::shared_ptr<GenericPart>;


} // namespace cie::gl


#endif