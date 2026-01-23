#pragma once

// --- Internal Includes ---
#include "poisson2D/MeshData.hpp"
#include "poisson2D/CellData.hpp"
#include "poisson2D/BoundaryData.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/LegendrePolynomial.hpp"

// ---- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"
#include "packages/trees/inc//ContiguousSpaceTree.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/commandline/inc/ArgParse.hpp"


namespace cie::fem {


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


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
void generateMesh(Ref<Mesh> rMesh,
                  Ref<const utils::ArgParse::Results> rArguments) {
    auto logBlock = utils::LoggerSingleton::get().newBlock("generate mesh");
    const unsigned nodesPerDirection = rArguments.get<std::size_t>("r");

    {
        // Define an ansatz space and its derivatives.
        // In this example, every cell will use the same ansatz space.
        Ansatz ansatzSpace;
        Ansatz::Derivative ansatzDerivative;

        std::array<Basis,polynomialOrder+1> basisFunctions;
        std::array<Basis::Derivative,polynomialOrder+1> basisDerivatives;

        {
            auto logBlock = utils::LoggerSingleton::get().newBlock("generate basis functions");

            // Construct basis functions.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                // Generate a 1D polynomial serving as one of the basis functions.
                maths::IntegratedLegendrePolynomial<Scalar> legendre(iBasis);
                std::array<Scalar,Basis::coefficientCount> polynomialCoefficients;
                CIE_CHECK(
                    legendre.coefficients().size() <= polynomialCoefficients.size(),
                    "basis function " << iBasis << " is expected to have at most "
                    << polynomialCoefficients.size() << " coefficients, but has "
                    << legendre.coefficients().size())
                std::copy_n(
                    legendre.coefficients().data(),
                    legendre.coefficients().size(),
                    polynomialCoefficients.data());
                std::fill_n(
                    polynomialCoefficients.data() + legendre.coefficients().size(),
                    legendre.coefficients().size() < polynomialCoefficients.size()
                        ? polynomialCoefficients.size() - legendre.coefficients().size()
                        : 0u,
                    0.0);

                // Construct a new polynomial view over the copied coefficients.
                // The array of coefficients is not stable, so these views must
                // be reassigned after basis generation is done.
                basisFunctions[iBasis] = Basis(polynomialCoefficients);
            } // for iBasis in range(polynomialOrder + 1)

            // Construct the derivatives of all basis functions.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                Ref<Basis> rBasis = basisFunctions[iBasis];
                basisDerivatives[iBasis] = rBasis.makeDerivative();
            } // for iBasis in range(polynomialOrder + 1)

            ansatzSpace = Ansatz(basisFunctions);
        }

        // Construct the derivatives of ansatz spaces.
        ansatzDerivative = Ansatz::Derivative(basisFunctions);

        rMesh.data() = MeshData(
            std::move(ansatzSpace),
            std::move(ansatzDerivative));
        }

    // Insert cells into the adjacency graph
    const Scalar edgeLength = 1.0 / (nodesPerDirection - 1);
    Size iBoundary = 0ul;

    for (Size iCellRow : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
        for (Size iCellColumn : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
            StaticArray<StaticArray<Scalar,Dimension>,2> transformed;
            OrientedAxes<Dimension> axes;

            // Define the cell's orientation in topological and physical space.
            if (iCellRow % 2) {
                axes[0] = "-x";
                transformed[0][0] = (iCellRow + 1.0) * edgeLength;
                transformed[1][0] = iCellRow * edgeLength;
            } else {
                axes[0] = "+x";
                transformed[0][0] = iCellRow * edgeLength;
                transformed[1][0] = (iCellRow + 1.0) * edgeLength;
            }

            if (iCellColumn % 2) {
                axes[1] = "-y";
                transformed[0][1] = (iCellColumn + 1.0) * edgeLength;
                transformed[1][1] = iCellColumn * edgeLength;
            } else {
                axes[1] = "+y";
                transformed[0][1] = iCellColumn * edgeLength;
                transformed[1][1] = (iCellColumn + 1.0) * edgeLength;
            }

            // Insert the cell into the adjacency graph (mesh) as a vertex
            const Size iCell = iCellRow * (nodesPerDirection - 1u) + iCellColumn;
            Mesh::Vertex::Data data (
                VertexID(iCell), // <= todo: remove duplicate id
                0u,   // <= All cells share the same ansatz space in this example.
                1.0,
                axes,
                SpatialTransform(transformed.begin(), transformed.end())
            );
            rMesh.insert(Mesh::Vertex(
                VertexID(iCell),
                {}, ///< edges of the adjacency graph are added automatically during edge insertion
                std::move(data)
            ));

            // Insert the current cell's connections to other cells already in the mesh.
            // The rule here is that cells with lower manhattan distance from the origin
            // are sources, while those with a higher norm are targets.
            if (iCellRow) {
                const Size iSourceCell = iCell - (nodesPerDirection - 1);
                const Size iTargetCell = iCell;
                BoundaryID sharedBoundary = iCellRow % 2 ? BoundaryID("+x") : BoundaryID("-x");
                rMesh.insert(Mesh::Edge(
                    EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    {sharedBoundary}
                ));
            } // if iCellRow

            if (iCellColumn) {
                const Size iSourceCell = iCell - 1ul;
                const Size iTargetCell = iCell;
                BoundaryID sharedBoundary = iCellColumn % 2 ? BoundaryID("+y") : BoundaryID("-y");
                rMesh.insert(Mesh::Edge(
                    EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    {sharedBoundary}
                ));
            } // if iCellColumn
        } // for iCellColumn in range(nodesPerDirection -1)
    } // for iCellRow in range(nodesPerDirection - 1)
}


} // namespace cie::fem
