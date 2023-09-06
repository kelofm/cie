#ifndef CIE_GL_BUFFER_VERTEX_3_HPP
#define CIE_GL_BUFFER_VERTEX_3_HPP

// --- Internal Includes ---
#include "packages/buffer/inc/Vertex2.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::gl {


class Vertex3 : public Vertex2
{
public:
    CIE_DEFINE_CLASS_POINTERS(Vertex3)

public:
    Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
            Size offset );

    Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
            Size offset,
            Vertex3::value_type x,
            Vertex3::value_type y,
            Vertex3::value_type z );

    virtual Size numberOfAttributes() const override;

    Vertex3::value_type& z();
    Vertex3::value_type z() const;

protected:
    /// Overloaded constructor for derived Vertex classes
    Vertex3(AttributeContainer::SharedPointer p_attributeContainer,
            Size offset,
            const Vertex3::offset_container& r_attributeOffsets );
};


} // namespace cie::gl

#endif
