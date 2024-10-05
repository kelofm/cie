#ifndef CIE_GL_GL_TRAITS_HPP
#define CIE_GL_GL_TRAITS_HPP

/// @defgroup ciegl Visualization Library (obsolete)

// --- External Includes ---
#include <glad/glad.h>

// --- Internal Includes ---
#include "packages/buffer/inc/AttributeContainer.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <vector>


namespace cie::gl {


struct GLTraits
{
    using attribute_type = AttributeContainer::value_type;

    using attribute_container = AttributeContainer;

    using index_type = GLuint;

    using index_container = std::vector<index_type>;

    using color_type = StaticArray<attribute_type,4>;
};


} // namespace cie::gl


#endif