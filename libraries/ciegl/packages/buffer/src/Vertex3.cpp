// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex3.hpp"


namespace cie::gl {


const Vertex3::offset_container vertex3OffsetContainer { 0 };


Vertex3::Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset,
                 const Vertex3::offset_container& r_attributeOffsets) :
    Vertex2(p_attributeContainer,
            offset,
            r_attributeOffsets)
{
    this->resizeAttributeContainer();
}


Vertex3::Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset) :
    Vertex3(p_attributeContainer,
            offset,
            vertex3OffsetContainer)
{
}


Vertex3::Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset,
                 Vertex3::value_type x,
                 Vertex3::value_type y,
                 Vertex3::value_type z) :
    Vertex3(p_attributeContainer, offset)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->x() = x;
    this->y() = y;
    this->z() = z;

    CIE_END_EXCEPTION_TRACING
}


Size Vertex3::numberOfAttributes() const
{
    return 3;
}


Vertex3::value_type& Vertex3::z()
{
    return this->at(2);
}


Vertex3::value_type Vertex3::z() const
{
    return this->at(2);
}


} // namespace cie::gl