// --- Internal Includes ---
#include "packages/buffer/inc/AttributeContainer.hpp"


namespace cie::gl {


AttributeContainer::AttributeContainer( const std::initializer_list<typename AttributeContainer::value_type>& r_list  ) :
    AttributeContainer::base_container( r_list )
{
}


} // namespace cie::gl