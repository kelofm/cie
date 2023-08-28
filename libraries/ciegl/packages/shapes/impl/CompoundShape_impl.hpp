#ifndef CIE_GL_SHAPES_COMPOUND_SHAPE_IMPL_HPP
#define CIE_GL_SHAPES_COMPOUND_SHAPE_IMPL_HPP

// --- Utility Includes ---
#include "packages/ranges/inc/TransformView.hpp"
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <numeric>


namespace cie::gl {


template <concepts::GLVertex VertexType>
CompoundShape<VertexType>::CompoundShape( typename CompoundShape<VertexType>::attribute_container::SharedPointer p_container ) :
    AbsGLShape( p_container ),
    _shapes()
{
}


template <concepts::GLVertex VertexType>
CompoundShape<VertexType>::CompoundShape( typename CompoundShape<VertexType>::shape_container&& r_shapes ) :
    AbsGLShape( r_shapes[0]->attributes() ),
    _shapes( std::forward<typename CompoundShape<VertexType>::shape_container>(r_shapes) )
{
    for ( const auto& rp_shape : this->_shapes )
        if ( rp_shape->attributes() != this->_shapes[0]->attributes() )
            CIE_THROW( Exception, "Attribute container mismatch!" )
}


template <concepts::GLVertex VertexType>
template <class ContainerType>
requires concepts::InterfaceContainer<ContainerType,GLShape<VertexType>>
CompoundShape<VertexType>::CompoundShape( ContainerType&& r_shapes ) :
    AbsGLShape( r_shapes[0]->attributes() ),
    _shapes()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_shapes.insert( this->_shapes.end(),
                          r_shapes.begin(),
                          r_shapes.end() );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
inline void
CompoundShape<VertexType>::setAttribute( Size attributeIndex,
                                         Size componentIndex,
                                         typename CompoundShape<VertexType>::attribute_type value )
{
    CIE_BEGIN_EXCEPTION_TRACING

    for ( auto& rp_shape : this->_shapes )
        rp_shape->setAttribute( attributeIndex,
                                componentIndex,
                                value );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
inline void
CompoundShape<VertexType>::updateShape()
{
}


template <concepts::GLVertex VertexType>
typename CompoundShape<VertexType>::index_container
CompoundShape<VertexType>::indices() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    typename CompoundShape<VertexType>::index_container indices, shapeIndices;
    Size indexOffset = 0;

    for ( const auto& rp_shape : this->_shapes )
    {
        shapeIndices = rp_shape->indices();
        utils::reserve( indices, indices.size() + shapeIndices.size() );

        auto offsetShapeIndices = utils::makeTransformView<Size>(
            shapeIndices,
            [indexOffset](auto index) { return index + indexOffset; }
        );

        indices.insert( indices.end(),
                        offsetShapeIndices.begin(),
                        offsetShapeIndices.end() );

        CIE_OUT_OF_RANGE_CHECK( !rp_shape->vertices().empty() )

        indexOffset += rp_shape->vertices().size();
    }

    return indices;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
inline const typename CompoundShape<VertexType>::shape_container&
CompoundShape<VertexType>::shapes() const
{
    return this->_shapes;
}


} // namespace cie::gl


#endif