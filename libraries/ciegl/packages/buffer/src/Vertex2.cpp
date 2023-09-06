// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex2.hpp"


namespace cie::gl {


const Vertex2::offset_container vertex2OffsetContainer { 0 };


Vertex2::Vertex2(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset,
                 const Vertex2::offset_container& r_attributeContainer ) :
    AbsVertex(p_attributeContainer,
              offset,
              r_attributeContainer )
{
    this->resizeAttributeContainer();
}


Vertex2::Vertex2(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset ) :
    Vertex2( p_attributeContainer,
             offset,
             vertex2OffsetContainer)
{
}


Vertex2::Vertex2(AttributeContainer::SharedPointer p_attributeContainer,
                 Size offset,
                 Vertex2::value_type x,
                 Vertex2::value_type y) :
    Vertex2(p_attributeContainer,
            offset)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->x() = x;
    this->y() = y;

    CIE_END_EXCEPTION_TRACING
}


Size Vertex2::numberOfAttributes() const
{
    return 2;
}


Vertex2::attribute_range Vertex2::position()
{
    return this->attribute(0);
}


Vertex2::attribute_const_range Vertex2::position() const
{
    return this->attribute(0);
}


Vertex2::value_type& Vertex2::x()
{
    return this->at(0);
}


Vertex2::value_type Vertex2::x() const
{
    return this->at(0);
}


Vertex2::value_type& Vertex2::y()
{
    return this->at(1);
}


Vertex2::value_type Vertex2::y() const
{
    return this->at(1);
}


} // namespace cie::gl