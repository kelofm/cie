#ifndef CIE_GL_SHAPES_AXES_HPP
#define CIE_GL_SHAPES_AXES_HPP

// --- Internal Includes ---
#include "packages/shapes/inc/CompoundShape.hpp"
#include "packages/buffer/inc/ColoredVertex3.hpp"
#include "packages/camera/inc/RigidBody.hpp"


namespace cie::gl {


class Axes :
    public CompoundShape<ColoredVertex3>,
    public RigidBody
{
public:
    Axes( Axes::attribute_container::SharedPointer p_attributes,
          Size numberOfAxes,
          const Axes::vector_type& r_position = {0.0, 0.0, 0.0},
          const Axes::vector_type& r_direction = {1.0, 0.0, 0.0},
          const Axes::vector_type& r_up = {0.0, 0.0, 1.0} );
};


} // namespace cie::gl


#endif