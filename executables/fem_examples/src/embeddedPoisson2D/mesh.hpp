#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/MeshData.hpp"
#include "embeddedPoisson2D/CellData.hpp"
#include "embeddedPoisson2D/BoundaryData.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"


namespace cie::fem {


/// @brief Mesh type.
/// @details The mesh is modeled by an adjacency graph whose vertices are cells,
///          and edges represent the boundaries between adjacent cells. The graph
///          is not directed, but the @ref BoundaryID of the graph edge is always
///          defined on the edge's source cell.
using Mesh = Graph<CellData,BoundaryData,MeshData>;


/// @brief Generate cells and boundaries for the example problem.
/// @details Mesh:
///          @code
///            [u(0,1)=2]                                     [u(1,1)=3]
///                     +---------+                 +---------+
///                     | (m-1)n+1|                 |    mn   |
///                     +---------+                 +---------+
///                       .
///                       .
///                       .
///                     +---------+
///                     |   n+1   |
///                     +---------+
///                     +---------+---------+       +---------+
///                     |    1    |    2    |  ...  |    n    |
///                     +---------+---------+       +---------+
///            [u(0,0)=0]                                     [u(1,0)=1]
///          @endcode
void generateMesh(
    Ref<Mesh> rMesh,
    std::span<const Scalar,2> meshBase,
    std::span<const Scalar,2> meshLengths,
    Ref<const utils::ArgParse::Results> rArguments);


} // namespace cie::fem
