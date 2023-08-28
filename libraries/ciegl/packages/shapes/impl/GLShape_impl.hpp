#ifndef CIE_GL_SHAPES_GL_SHAPE_IMPL_HPP
#define CIE_GL_SHAPES_GL_SHAPE_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <concepts::GLVertex VertexType>
GLShape<VertexType>::GLShape( typename GLShape<VertexType>::attribute_container::SharedPointer p_attributes ) :
    AbsGLShape( p_attributes ),
    _vertices()
{
}


template <concepts::GLVertex VertexType>
void GLShape<VertexType>::setAttribute( Size attributeIndex,
                                        Size componentIndex,
                                        typename GLShape<VertexType>::attribute_type value )
{
    CIE_BEGIN_EXCEPTION_TRACING

    for ( auto& rp_vertex : this->_vertices )
        rp_vertex->at( attributeIndex, componentIndex ) = value;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
inline const typename GLShape<VertexType>::vertex_container&
GLShape<VertexType>::vertices() const
{
    return this->_vertices;
}


} // namespace cie::gl


#endif