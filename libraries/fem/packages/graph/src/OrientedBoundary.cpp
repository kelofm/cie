// --- FEM Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"

// --- STL Includes ---
#include <iostream> // ostream


namespace cie::fem {


template <unsigned Dimension>
std::ostream& operator<<(std::ostream& rStream, OrientedBoundary<Dimension> boundary)
{
    for (auto component : boundary) rStream << component;
    rStream << ':' << boundary.id();
    return rStream;
}


#define CIE_DEFINE_ORIENTED_BOUNDARY_SERIALIZATION(Dimension) \
    template std::ostream& operator<<(std::ostream&,OrientedBoundary<Dimension>)

CIE_DEFINE_ORIENTED_BOUNDARY_SERIALIZATION(1);
CIE_DEFINE_ORIENTED_BOUNDARY_SERIALIZATION(2);
CIE_DEFINE_ORIENTED_BOUNDARY_SERIALIZATION(3);


} // namespace cie::fem
