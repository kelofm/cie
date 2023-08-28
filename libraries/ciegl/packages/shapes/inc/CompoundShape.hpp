#ifndef CIE_GL_SHAPES_COMPOUND_SHAPE_HPP
#define CIE_GL_SHAPES_COMPOUND_SHAPE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/shapes/inc/GLShape.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>


namespace cie::gl {


/// Container of shapes with an AbsGLShape interface
template <concepts::GLVertex VertexType>
class CompoundShape : public AbsGLShape
{
public:
    using vertex_type     = VertexType;
    using shape_type      = GLShape<vertex_type>;
    using shape_ptr       = std::shared_ptr<shape_type>;
    using shape_container = std::vector<shape_ptr>;

    using AbsGLShape::attribute_type;
    using AbsGLShape::index_container;

public:
    CompoundShape( attribute_container::SharedPointer p_container );

    /// @note will segfault on an empty container!
    CompoundShape( shape_container&& r_shapes );

    /// @note will segfault on an empty container!
    template <class ContainerType>
    requires concepts::InterfaceContainer<ContainerType,GLShape<VertexType>>
    CompoundShape( ContainerType&& r_shapes );

    virtual void setAttribute( Size attributeIndex,
                               Size componentIndex,
                               attribute_type value ) override;

    virtual void updateShape() override;

    virtual index_container indices() const override;

    const shape_container& shapes() const;

protected:
    shape_container _shapes;
};


} // namespace cie::gl

#include "packages/shapes/impl/CompoundShape_impl.hpp"

#endif