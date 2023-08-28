#ifndef CIE_GL_SHAPES_ARROW_IMPL_HPP
#define CIE_GL_SHAPES_ARROW_IMPL_HPP

#define _USE_MATH_INCLUDES

// --- LinAlg Includes ---
#include "packages/overloads/inc/vectoroperators.hpp"
#include "packages/types/inc/typeoperations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/resize.hpp"

// --- STL Includes ---
#include <cmath>
#include <math.h>


namespace cie::gl {


template <concepts::GLVertex VertexType>
Arrow<VertexType>::Arrow(typename Arrow<VertexType>::attribute_container::SharedPointer p_attributes,
                         const typename Arrow<VertexType>::Point& r_source,
                         const typename Arrow<VertexType>::Point& r_sink,
                         typename Arrow<VertexType>::attribute_type headRatio,
                         typename Arrow<VertexType>::attribute_type baseRadius,
                         typename Arrow<VertexType>::attribute_type headRadius,
                         const Size resolution) :
    GLShape<VertexType>( p_attributes ),
    _source( r_source ),
    _sink( r_sink ),
    _headRatio( headRatio ),
    _baseRadius( baseRadius ),
    _headRadius( headRadius ),
    _resolution( resolution )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check arguments
    CIE_CHECK_POINTER( p_attributes )
    CIE_CHECK( 1 < resolution, "Resolution of a circle must be greater than 1" )

    // Allocate vertices
    Size numberOfVertices = 3 * resolution + 1 + 1;

    // Source vertex ( VertexType::numberOfAttributes needs an existing instance )
    auto p_source = makeVertex<VertexType>( p_attributes );
    this->_vertices.push_back( p_source );

    utils::reserve( *p_attributes,
                    p_attributes->size() + (numberOfVertices-1) * p_source->numberOfAttributes() );

    for ( Size vertexIndex=1; vertexIndex<numberOfVertices; ++vertexIndex )
        this->_vertices.push_back( makeVertex<VertexType>( p_attributes ) );

    // Compute shape
    this->updateShape();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
void
Arrow<VertexType>::updateShape()
{
    CIE_BEGIN_EXCEPTION_TRACING

    using attribute_type = typename Arrow<VertexType>::attribute_type;
    using Point = typename Arrow<VertexType>::Point;

    // Pull overloads from linalg
    using ::operator+;
    using ::operator-;
    using ::operator*;

    auto tangent = this->_sink - this->_source;

    attribute_type tolerance = 1e-8;

    // Check shape parameters
    CIE_DEBUG_CHECK(
        tolerance < linalg::norm( tangent ),
        "Degenerate arrow" )

    CIE_CHECK(
        0.0 < this->_headRatio && this->_headRatio < 1.0,
        "Invalid arrow head ratio" )

    CIE_CHECK(
        0.0 < this->_baseRadius,
        "Invalid arrow base radius" )

    CIE_CHECK(
        0.0 < this->_headRadius,
        "Invalid arrow head radius" )

    // Local definition of a cross product
    // TODO: implement generaé cross product somewhere else,
    // probably in linalg
    auto cross = []( const Point& r_lhs, const Point& r_rhs ) -> Point
    {
        return {
            r_lhs[1] * r_rhs[2] - r_rhs[1] * r_lhs[2],
            r_rhs[0] * r_lhs[2] - r_lhs[0] * r_rhs[2],
            r_lhs[0] * r_rhs[1] - r_rhs[0] * r_lhs[1]
        };
    };

    // Get basis on the plane normal to the arrow
    Point normal;

    if ( std::abs(tangent[1]) > tolerance )
        normal = cross( {1.0, 0.0, 0.0}, tangent );
    else
        normal = cross( {0.0, 1.0, 0.0}, tangent );

    linalg::normalize( normal, &linalg::norm<Point>, tolerance );

    Point binormal = cross( tangent, normal );

    linalg::normalize( binormal, &linalg::norm<Point>, tolerance );

    // <-- normal and binormal are unit vectors, the tangent is not

    // Functor for setting the position of a vertex
    auto setVertexPosition = []( typename Arrow<VertexType>::vertex_ptr p_vertex, const Point& r_position ) -> void
    {
        CIE_BEGIN_EXCEPTION_TRACING

        auto vertexPosition = p_vertex->position();
        for ( Size dim=0; dim<vertexPosition.size(); ++dim )
            vertexPosition[dim] = r_position[dim];

        CIE_END_EXCEPTION_TRACING
    };

    // Set position of the source vertex
    setVertexPosition(
        this->_vertices[0],
        this->_source
    );

    // Set position of the vertices on the circles
    attribute_type angleIncrement = 2.0 * M_PI / attribute_type(this->_resolution);

    attribute_type baseRatio = 1.0 - this->_headRatio;

    for ( Size vertexOnCircleIndex=0; vertexOnCircleIndex<this->_resolution; ++vertexOnCircleIndex )
    {
        attribute_type angle               = vertexOnCircleIndex * angleIncrement;
        attribute_type normalCoefficient   = std::sin( angle );
        attribute_type binormalCoefficient = std::cos( angle );
        Point unitOnPlane         = normalCoefficient*normal + binormalCoefficient*binormal;

        // Base - source end
        setVertexPosition(
            this->_vertices[vertexOnCircleIndex + 1],
            this->_source + this->_baseRadius * unitOnPlane
        );

        // Base - head end
        setVertexPosition(
            this->_vertices[this->_resolution + vertexOnCircleIndex + 1],
            this->_source + baseRatio * tangent + this->_baseRadius * unitOnPlane
        );

        // Head
        setVertexPosition(
            this->_vertices[2 * this->_resolution + vertexOnCircleIndex + 1],
            this->_source + baseRatio * tangent + this->_headRadius * unitOnPlane
        );
    }

    // Set tip position
    setVertexPosition(
        this->_vertices[3 * this->_resolution + 1],
        this->_sink
    );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::GLVertex VertexType>
typename Arrow<VertexType>::index_container
Arrow<VertexType>::indices() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    typename Arrow<VertexType>::index_container indices;

    Size numberOfTriangles = 6 * this->_resolution;

    utils::reserve( indices, 3 * numberOfTriangles );

    // Triangles on the base - source end
    for ( Size localIndex=0; localIndex<this->_resolution; ++localIndex )
    {
        indices.push_back( 0 ); // <-- source
        indices.push_back( 1 + localIndex );
        indices.push_back( 1 + (localIndex+1) % this->_resolution );
    }

    // Triangles on the base hull
    for ( Size localIndex=0; localIndex<this->_resolution; ++localIndex )
    {
        Size offset = 1;

        indices.push_back( offset + localIndex );
        indices.push_back( offset + this->_resolution + localIndex );
        indices.push_back( offset + this->_resolution + (localIndex+1) % this->_resolution );

        indices.push_back( offset + localIndex );
        indices.push_back( offset + this->_resolution + (localIndex+1) % this->_resolution );
        indices.push_back( offset + (localIndex+1) % this->_resolution );
    }

    // Triangles on the base-head interface
    for ( Size localIndex=0; localIndex<this->_resolution; ++localIndex )
    {
        Size offset = 1 + this->_resolution;

        indices.push_back( offset + localIndex );
        indices.push_back( offset + this->_resolution + localIndex );
        indices.push_back( offset + this->_resolution + (localIndex+1) % this->_resolution );

        indices.push_back( offset + localIndex );
        indices.push_back( offset + this->_resolution + (localIndex+1) % this->_resolution );
        indices.push_back( offset + (localIndex+1) % this->_resolution );
    }

    // Triangles on the head
    for ( Size localIndex=0; localIndex<this->_resolution; ++localIndex )
    {
        Size offset = 1 + 2 * this->_resolution;
        Size tipIndex = this->_vertices.size() - 1;

        indices.push_back( offset + localIndex );
        indices.push_back( tipIndex );
        indices.push_back( offset + (localIndex+1) % this->_resolution );
    }

    return indices;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif