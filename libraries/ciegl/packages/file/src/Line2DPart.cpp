// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/file/inc/Line2DPart.hpp"


namespace cie::gl {


void Line2DPart::load(std::istream&)
{
    CIE_THROW( NotImplementedException, "No loading routine is specified for Line2DPart" )
}


Size Line2DPart::dimension() const
{
    return 2;
}


Size Line2DPart::vertexAttributeSize() const
{
    return 2;
}


Size Line2DPart::primitiveVertexSize() const
{
    return 2;
}


void Line2DPart::writeSTL(const std::filesystem::path&) const
{
    CIE_THROW( Exception, "Cannot triangulate Line2DPart" )
}


} // namespace cie::gl