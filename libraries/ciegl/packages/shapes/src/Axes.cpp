// --- LinAlg Includes ---
#include "packages/overloads/inc/vectoroperators.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/shapes/inc/Axes.hpp"
#include "packages/shapes/inc/Arrow.hpp"


namespace cie::gl {


Axes::Axes(Axes::attribute_container::SharedPointer p_attributes,
           Size numberOfAxes,
           const Axes::vector_type& r_position,
           const Axes::vector_type& r_direction,
           const Axes::vector_type& r_up) :
    CompoundShape<Axes::vertex_type>( p_attributes ),
    RigidBody( r_position,
               r_direction,
               r_up )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(
        1 <= numberOfAxes && numberOfAxes <= 3,
        "Axes are implemented only for dimensions 1, 2, and 3"
    )

    // Pull in vector addition operator from linalg
    using ::operator+;

    // Convenience typedefs
    using ArrowType  = Arrow<Axes::vertex_type>;
    using Point = ArrowType::Point;

    // Convert position to types expected as arguments
    Point position {
        Point::value_type( r_position[0] ),
        Point::value_type( r_position[1] ),
        Point::value_type( r_position[2] )
    };

    this->_shapes.reserve( numberOfAxes );

    // Add x axis (red)
    this->_shapes.emplace_back( new ArrowType(
        p_attributes,
        position,
        position + Point {1.0, 0.0, 0.0}
    ) );

    this->_shapes.back()->setAttribute( 1, 1, 0.0 );    // set green attribute
    this->_shapes.back()->setAttribute( 1, 2, 0.0 );    // set blue attribute

    if ( 1 < numberOfAxes )
    {
        // Add y axis (green)
        this->_shapes.emplace_back( new ArrowType(
            p_attributes,
            position,
            position + Point {0.0, 1.0, 0.0}
        ) );

        this->_shapes.back()->setAttribute( 1, 0, 0.0 );    // set red attribute
        this->_shapes.back()->setAttribute( 1, 2, 0.0 );    // set blue attribute

        if ( 2 < numberOfAxes )
        {
            // Add z axis (blue)
            this->_shapes.emplace_back( new ArrowType(
                p_attributes,
                position,
                position + Point {0.0, 0.0, 1.0}
            ) );

            this->_shapes.back()->setAttribute( 1, 0, 0.0 );    // set red attribute
            this->_shapes.back()->setAttribute( 1, 1, 0.0 );    // set green attribute
        }
    }

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl