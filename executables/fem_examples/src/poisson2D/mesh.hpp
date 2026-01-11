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
                  Ref<const utils::ArgParse::Results> rArguments)
{
    auto logBlock = utils::LoggerSingleton::get().newBlock("generate mesh");

    const unsigned polynomialOrder = rArguments.get<std::size_t>("p");
    const unsigned nodesPerDirection = rArguments.get<std::size_t>("r");

    {
        // Define an ansatz space and its derivatives.
        // In this example, every cell will use the same ansatz space.
        DynamicArray<Scalar> polynomialCoefficients;
        DynamicArray<Basis> basisFunctions;
        DynamicArray<Basis::Derivative> basisDerivatives;
        DynamicArray<Ansatz> ansatzSpaces;
        DynamicArray<AnsatzDerivative> ansatzDerivatives;

        {
            auto logBlock = utils::LoggerSingleton::get().newBlock("generate basis functions");

            // Construct basis functions.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                // Generate a 1D polynomial serving as one of the basis functions.
                maths::IntegratedLegendrePolynomial<Scalar> legendre(iBasis);

                // Copy its coefficients to the array of all polynomial coefficients.
                const std::size_t iCoefficientBegin = polynomialCoefficients.size();
                polynomialCoefficients.resize(polynomialCoefficients.size() + legendre.coefficients().size());
                std::copy_n(
                    legendre.coefficients().data(),
                    legendre.coefficients().size(),
                    polynomialCoefficients.data() + iCoefficientBegin);

                // Construct a new polynomial view over the copied coefficients.
                // The array of coefficients is not stable, so these views must
                // be reassigned after basis generation is done.
                basisFunctions.emplace_back(std::span<const Scalar>(
                    polynomialCoefficients.data() + iCoefficientBegin,
                    legendre.coefficients().size()));
            } // for iBasis in range(polynomialOrder + 1)

            if (basisFunctions.size() < 2 * polynomialCoefficients.size()) {
                polynomialCoefficients.resize(2 * polynomialCoefficients.size() - basisFunctions.size());
            }

            // Reassign the coefficient ranges of all basis functions,
            // now that the list of polynomial coefficients is stable.
            std::size_t iBegin = 0ul;
            for (auto& rFunction : basisFunctions) {
                const std::size_t coefficientCount = rFunction.coefficients().size();
                rFunction = Basis({
                    polynomialCoefficients.data() + iBegin,
                    coefficientCount});
                iBegin += coefficientCount;
            } // for rFunction : basisFunctions

            // Construct the derivatives of all basis functions.
            // Each derivative will need 1 less coefficient.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                Ref<Basis> rBasis = basisFunctions[iBasis];
                const std::size_t coefficientCount = rBasis.coefficients().empty()
                    ? 0ul
                    : rBasis.coefficients().size() - 1;

                // Extend basis derivatives.
                basisDerivatives.emplace_back(rBasis.makeDerivative({
                    polynomialCoefficients.data() + iBegin,
                    coefficientCount}));

                iBegin += coefficientCount;
            } // for iBasis in range(polynomialOrder + 1)

            ansatzSpaces.emplace_back(basisFunctions);
        }

        // Construct the derivatives of ansatz spaces.
        ansatzDerivatives.emplace_back(
            basisFunctions,
            basisDerivatives);

        // Define integration orders.
        StaticArray<unsigned,1> integrationOrders;
        integrationOrders.front() = rArguments.get<std::size_t>("i");

        rMesh.data() = MeshData(
            integrationOrders,
            std::move(polynomialCoefficients),
            std::move(basisFunctions),
            std::move(basisDerivatives),
            std::move(ansatzSpaces),
            std::move(ansatzDerivatives));
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
