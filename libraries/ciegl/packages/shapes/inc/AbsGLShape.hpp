#ifndef CIE_GL_SHAPES_ABS_GL_SHAPE_HPP
#define CIE_GL_SHAPES_ABS_GL_SHAPE_HPP

// --- Internal Includes ---
#include "packages/traits/inc/GLTraits.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>


namespace cie::gl {


class AbsGLShape : public GLTraits
{
public:
    using GLTraits::index_type;
    using GLTraits::index_container;
    using GLTraits::attribute_type;
    using GLTraits::attribute_container;

public:
    AbsGLShape(attribute_container::SharedPointer _p_attributes);

    AbsGLShape( AbsGLShape&& r_rhs ) = default;

    virtual void setAttribute(Size attributeIndex,
                              Size componentIndex,
                              attribute_type value) = 0;

    /// Modify vertex attributes
    virtual void updateShape() = 0;

    /// Get vertex connectivity
    virtual index_container indices() const = 0;

    /// Get pointer to attribute container
    const attribute_container::SharedPointer attributes() const;

    /// Get pointer to attribute container
    attribute_container::SharedPointer attributes();

protected:
    AbsGLShape() = delete;
    AbsGLShape( const AbsGLShape& r_rhs ) = delete;
    AbsGLShape& operator=( const AbsGLShape& r_rhs ) = delete;

protected:
    attribute_container::SharedPointer _p_attributes;
};


using GLShapePtr = std::shared_ptr<AbsGLShape>;


} // namespace cie::gl

#endif