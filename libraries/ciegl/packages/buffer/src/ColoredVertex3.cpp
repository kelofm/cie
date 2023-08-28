// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/ColoredVertex3.hpp"

// --- STL Includes ---
#include <algorithm>
#include <numeric>


namespace cie::gl {


const ColoredVertex3::offset_container coloredVertex3OffsetContainer = { 0, 3 };


ColoredVertex3::ColoredVertex3(AttributeContainer::SharedPointer p_attributeContainer,
                               Size offset,
                               const ColoredVertex3::Point& r_point,
                               const GLTraits::color_type& r_color) :
    Vertex3( p_attributeContainer,
             offset,
             coloredVertex3OffsetContainer )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->resizeAttributeContainer();

    std::copy( r_point.begin(),
               r_point.end(),
               this->attribute(0).begin() );

    std::copy( r_color.begin(),
               r_color.end(),
               this->attribute(1).begin() );

    CIE_END_EXCEPTION_TRACING
}


ColoredVertex3::ColoredVertex3( AttributeContainer::SharedPointer p_attributeContainer,
                                Size offset ) :
    ColoredVertex3( p_attributeContainer,
                    offset,
                    {0.0, 0.0, 0.0},
                    {1.0, 1.0, 1.0, 1.0} )
{
}


const Size ColoredVertex3::numberOfAttributes() const
{
    return 7;
}


ColoredVertex3::attribute_range ColoredVertex3::color()
{
    CIE_BEGIN_EXCEPTION_TRACING

    return this->attribute(1);

    CIE_END_EXCEPTION_TRACING
}


ColoredVertex3::attribute_const_range ColoredVertex3::color() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return this->attribute(1);

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl
