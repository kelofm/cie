#ifndef CIE_GL_BUFFER_VERTEX_2_HPP
#define CIE_GL_BUFFER_VERTEX_2_HPP

// --- Internal Includes ---
#include "packages/buffer/inc/AbsVertex.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::gl {


class Vertex2 : public AbsVertex
{
public:
    CIE_DEFINE_CLASS_POINTERS(Vertex2)

public:
    Vertex2(AttributeContainer::SharedPointer p_attributeContainer,
            Size offset );

    Vertex2( AttributeContainer::SharedPointer p_attributeContainer,
             Size offset,
             Vertex2::value_type x,
             Vertex2::value_type y );

    virtual const Size numberOfAttributes() const override;

    virtual Vertex2::attribute_range position();
    virtual Vertex2::attribute_const_range position() const;

    Vertex2::value_type& x();
    const Vertex2::value_type x() const;

    Vertex2::value_type& y();
    const Vertex2::value_type y() const;

protected:
    /// Overloaded constructor for derived Vertex classes
    Vertex2(AttributeContainer::SharedPointer p_attributeContainer,
            Size offset,
            const Vertex2::offset_container& r_attributeOffsets);
};


} // namespace cie::gl

#endif