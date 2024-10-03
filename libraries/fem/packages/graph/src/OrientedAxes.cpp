// --- FEM Includes ---
#include "packages/graph/inc/OrientedAxes.hpp"

// --- STL Includes ---
#include <iostream> // ostream


namespace cie::fem {


template <unsigned Dimension>
std::ostream& operator<<(std::ostream& rStream, OrientedAxes<Dimension> boundary)
{
    for (auto component : boundary) rStream << component;
    return rStream;
}


#define CIE_DEFINE_ORIENTED_AXES_SERIALIZATION(Dimension) \
    template std::ostream& operator<<(std::ostream&,OrientedAxes<Dimension>)

CIE_DEFINE_ORIENTED_AXES_SERIALIZATION(1);
CIE_DEFINE_ORIENTED_AXES_SERIALIZATION(2);
CIE_DEFINE_ORIENTED_AXES_SERIALIZATION(3);


} // namespace cie::fem
