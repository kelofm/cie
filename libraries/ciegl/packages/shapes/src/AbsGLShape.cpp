// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/shapes/inc/AbsGLShape.hpp"


namespace cie::gl {



AbsGLShape::AbsGLShape( AbsGLShape::attribute_container::SharedPointer p_attributes ) :
    _p_attributes( p_attributes )
{
}


const AbsGLShape::attribute_container::SharedPointer AbsGLShape::attributes() const
{
    return this->_p_attributes;
}


AbsGLShape::attribute_container::SharedPointer AbsGLShape::attributes()
{
    return this->_p_attributes;
}


} // namespace cie::gl