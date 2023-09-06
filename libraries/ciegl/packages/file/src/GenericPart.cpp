// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/file/inc/GenericPart.hpp"


namespace cie::gl {


GenericPart::GenericPart(const Size dimension,
                         const Size vertexAttributeSize,
                         const Size primitiveVertexSize) :
    _dimension(dimension),
    _vertexAttributeSize(vertexAttributeSize),
    _primitiveVertexSize(primitiveVertexSize)
{
}


GenericPart::attribute_container::SharedPointer GenericPart::attributes()
{
    return this->_p_attributes;
}


GenericPart::index_container& GenericPart::indices()
{
    return this->_indices;
}


Size GenericPart::dimension() const
{
    return this->_dimension;
}


Size GenericPart::vertexAttributeSize() const
{
    return this->_vertexAttributeSize;
}


Size GenericPart::primitiveVertexSize() const
{
    return this->_primitiveVertexSize;
}


void GenericPart::writeSTL(const std::filesystem::path&) const
{
    CIE_THROW( NotImplementedException, "Cannot write GenericPart to STL" )
}


} // namespace cie::gl