// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"


namespace cie::fem {


Ref<std::ostream> operator<<(Ref<std::ostream> rStream, BoundaryID id)
{
    rStream << (id.getDirection() ? '+' : '-') << id.getDimension();
    return rStream;
}


} // namespace cie::fem
