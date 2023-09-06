#ifndef CIE_GL_BUFFER_COLORED_VERTEX_3_HPP
#define CIE_GL_BUFFER_COLORED_VERTEX_3_HPP

// --- Utility Includes ---
#include "packages/ranges/inc/TransformView.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex3.hpp"
#include "packages/traits/inc/GLTraits.hpp"


namespace cie::gl {


class ColoredVertex3 : public Vertex3
{
public:
    using Point = StaticArray<ColoredVertex3::value_type,3>;
    using view_type  = utils::TransformView<AttributeContainer,AttributeContainer::value_type>;
    using color_type = GLTraits::color_type;

public:
    ColoredVertex3(AttributeContainer::SharedPointer p_attributeContainer,
                   Size offset,
                   const Point& r_point,
                   const GLTraits::color_type& r_color);

    ColoredVertex3(AttributeContainer::SharedPointer p_attributeContainer,
                   Size offset);

    virtual Size numberOfAttributes() const override;

    ColoredVertex3::attribute_range color();

    ColoredVertex3::attribute_const_range color() const;

protected:
    ColoredVertex3() = delete;

    /// Constructor for derived classes
    ColoredVertex3(AttributeContainer::SharedPointer p_attributeContainer,
                   Size offset,
                   const ColoredVertex3::offset_container& r_attributeOffset);
};


} // namespace cie::gl


#endif
