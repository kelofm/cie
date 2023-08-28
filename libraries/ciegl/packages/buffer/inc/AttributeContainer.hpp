#ifndef CIE_GL_BUFFER_ATTRIBUTE_CONTAINER_HPP
#define CIE_GL_BUFFER_ATTRIBUTE_CONTAINER_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>
#include <concepts>


namespace cie::gl {


class AbsVertex;
using VertexPtr = std::shared_ptr<AbsVertex>;


class AttributeContainer :
    public DynamicArray<GLfloat>,
    public std::enable_shared_from_this<AttributeContainer>
{
public:
    using base_container = DynamicArray<GLfloat>;

    CIE_DEFINE_CLASS_POINTERS(AttributeContainer)

public:
    AttributeContainer( const std::initializer_list<typename AttributeContainer::value_type>& r_list );

    template <class ...Args>
    AttributeContainer( Args&&... args );
};


} // namespace cie::gl

#include "packages/buffer/impl/AttributeContainer_impl.hpp"

#endif