#ifndef CIE_GL_SHAPES_GL_SHAPE_HPP
#define CIE_GL_SHAPES_GL_SHAPE_HPP

// --- Internal Includes ---
#include "packages/shapes/inc/AbsGLShape.hpp"
#include "packages/buffer/inc/AbsVertex.hpp"


namespace cie::gl {


template <concepts::GLVertex VertexType>
class GLShape : public AbsGLShape
{
public:
    using vertex_type      = VertexType;

    using vertex_ptr       = std::shared_ptr<vertex_type>;

    using vertex_container = std::vector<vertex_ptr>;

    using AbsGLShape::attribute_type;

    using AbsGLShape::attribute_container;

    using AbsGLShape::index_type;

    using AbsGLShape::index_container;

public:
    GLShape(attribute_container::SharedPointer p_attributes);

    GLShape( GLShape<VertexType>&& r_rhs ) = default;

    virtual void setAttribute( Size attributeIndex,
                               Size componentIndex,
                               attribute_type value ) override;

    const vertex_container& vertices() const;

protected:
    GLShape() = delete;

    GLShape( const GLShape<VertexType>& r_rhs ) = delete;

    GLShape<VertexType>& operator=( const GLShape<VertexType>& r_rhs ) = delete;

protected:
    vertex_container _vertices;
};


} // namespace cie::gl

#include "packages/shapes/impl/GLShape_impl.hpp"

#endif