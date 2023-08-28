#ifndef CIE_GL_SHAPES_ARROW_HPP
#define CIE_GL_SHAPES_ARROW_HPP

// --- Internal Includes ---
#include "packages/shapes/inc/GLShape.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <array>


namespace cie::gl {


template <concepts::GLVertex VertexType>
class Arrow : public GLShape<VertexType>
{
public:
    using typename GLShape<VertexType>::attribute_type;
    using typename GLShape<VertexType>::index_container;

    using Point = StaticArray<attribute_type,3>;

public:
    /**
     * @param p_attributes pointer to the attribute container in which the vertex attributes are stored
     * @param r_source source point
     * @param r_sink sin point
     * @param headRatio ratio of the head length to the total length
     * @param baseRadius radius of the base
     * @param headRadius largest radius of the arrowhead
     * @param resolution number of vertices on the base circles
     */
    Arrow(typename Arrow<VertexType>::attribute_container::SharedPointer p_attributes,
          const Point& r_source,
          const Point& r_sink,
          attribute_type headRatio = 0.4,
          attribute_type baseRadius = 0.1,
          attribute_type headRadius = 0.2,
          const Size resolution = 4);

    Arrow( Arrow<VertexType>&& r_rhs ) = default;

    virtual void updateShape() override;

    virtual index_container indices() const override;

protected:
    Point    _source;
    Point    _sink;
    attribute_type _headRatio;
    attribute_type _baseRadius;
    attribute_type _headRadius;
    const Size    _resolution;
};


} // namespace cie::gl

#include "packages/shapes/impl/Arrow_impl.hpp"

#endif