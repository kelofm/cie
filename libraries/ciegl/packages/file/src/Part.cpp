// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/file/inc/Part.hpp"

// --- STL Includes ---
#include <fstream>


namespace cie::gl {


Part::Part() :
    _p_attributes( new Part::attribute_container ),
    _indices()
{
}


Size Part::attributeByteCount() const
{
    return this->_p_attributes->size() * sizeof( Part::attribute_type );
}


Size Part::indexByteCount() const
{
    return this->_indices.size() * sizeof( Part::index_type );
}


Size Part::numberOfVertices() const
{
    return this->_p_attributes->size() / this->vertexAttributeSize();
}


Size Part::numberOfPrimitives() const
{
    return this->_indices.size() / this->primitiveVertexSize();
}


void Part::load(std::istream&)
{
    CIE_THROW( AbstractCallException, "Call to pure virtual Part::load" )
}


const Part::attribute_container::SharedPointer& Part::attributes() const
{
    return this->_p_attributes;
}


const Part::index_container& Part::indices() const
{
    return this->_indices;
}


Size Part::primitiveAttributeSize() const
{
    return this->vertexAttributeSize() * this->primitiveVertexSize();
}


Size Part::primitiveByteSize() const
{
    return this->primitiveAttributeSize() * sizeof( Part::attribute_type );
}


} // namespace cie::gl