#ifndef CIE_GL_BUFFER_ATTRIBUTE_CONTAINER_IMPL_HPP
#define CIE_GL_BUFFER_ATTRIBUTE_CONTAINER_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/resize.hpp"


namespace cie::gl {


template <class ...Args>
AttributeContainer::AttributeContainer( Args&&... args ) :
    AttributeContainer::base_container( std::forward<Args>(args)... )
{
}


} // namespace cie::gl

#endif