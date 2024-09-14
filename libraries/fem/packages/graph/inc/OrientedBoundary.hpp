#ifndef CIE_FEM_ORIENTED_BOUNDARY_HPP
#define CIE_FEM_ORIENTED_BOUNDARY_HPP

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"

// --- STL Includes ---
#include <bitset>


namespace cie::fem {


namespace detail {
template <unsigned Dimension>
class CubeVertexSet
{
private:
    constexpr static unsigned getIntegerBitSize()
    {
        unsigned iMax = intPow(2u, Dimension) - 1;
        unsigned size = 0u;
        while (iMax) {
            iMax >>= 1;
            ++size;
        }
        return size;
    }

public:
    static constexpr unsigned integerBitSize = getIntegerBitSize();

    static constexpr unsigned cornerSetBitSize = integerBitSize * intPow(2u, Dimension);
}; // struct CubeVertexSet
} // namespace detail


/** @brief Class defining a $d-1$-diensional cube spanning $[-1, 1]^{d-1}$ in any axis-aligned rotation.

 */
template <unsigned Dimension>
class OrientedBoundary
{
public:

private:
    BoundaryID _id;

    using Vertices = std::bitset<std::conditional_t<
        Dimension == 0u,
        detail::CubeVertexSet<0u>,
        detail::CubeVertexSet<Dimension - 1u>
    >::cornerSetBitSize>;
    Vertices _vertices;
}; // class OrientedBoundary


} // namespace cie::fem


#endif
