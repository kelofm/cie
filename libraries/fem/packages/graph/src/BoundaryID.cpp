// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"


namespace cie::fem {


Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, BoundaryID id)
{
    r_stream << (id.getDirection() ? '+' : '-') << id.getDimension();
    return r_stream;
}


} // namespace cie::fem
