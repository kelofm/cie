// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/MatrixMarket.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/utilities/inc/Graphviz.hpp"
#include "packages/maths/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/maths/inc/TransformedIntegrand.hpp"

// --- STL Includes ---
#include <ranges> // ranges::iota


namespace cie::fem {


using Scalar = double;


template <class TAnsatzSpace>
struct CellData
{
    using SpatialTransform = maths::ScaleTranslateTransform<Scalar,2u>;

    Scalar                                              diffusivity;
    OrientedAxes<TAnsatzSpace::Dimension>               axes;
    SpatialTransform                                    spatialTransform;
    std::shared_ptr<TAnsatzSpace>                       pAnsatzSpace;
    std::shared_ptr<typename TAnsatzSpace::Derivative>  pAnsatzDerivatives;
}; // struct CellData


/** @brief 2D system test.
 *  @details Mesh:
 *           @code
 *                  [u(0,1)=2]          [u(1,1)=3]
 *                      +---+           +---+
 *                      | 0 |           | 0 |
 *                      +---+           +---+
 *                        .
 *                        .
 *                        .
 *                      +---+
 *                      |n+1|
 *                      +---+
 *                      +---+---+       +---+
 *                      | 0 | 1 |  ...  | n |
 *                      +---+---+       +---+
 *                  [u(0,0)=0]          [u(1,0)=1]
 *           @endcode
 */
CIE_TEST_CASE("2D", "[systemTests]")
{
    CIE_TEST_CASE_INIT("2D")
    constexpr unsigned Dimension = 2;
    constexpr Size nodesPerDirection = 11;

    // Construct a 2D ansatz space.
    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;
    auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
        ,Basis({ 1.0,  0.0, -1.0})
    });
    auto pAnsatzDerivatives = std::make_shared<Ansatz::Derivative>(pAnsatzSpace->makeDerivative());

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent elements, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const StaticArray<Scalar,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    const auto ansatzMap = makeAnsatzMap(*pAnsatzSpace,
                                         samples,
                                         utils::Comparison<Scalar>(1e-8, 1e-6));

    // Construct mesh
    // The mesh is modeled by an adjacency graph whose vertices are cells,
    // and edges represent the boundaries between adjacent cells. The graph
    // is not directed, but the BoundaryID of the graph edge is always
    // defined on the edge's source cell.
    using Mesh = Graph<
        CellData<Ansatz>,
        std::pair<BoundaryID,BoundaryID>
    >;
    Mesh mesh;

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
            mesh.insert(Mesh::Vertex(
                Mesh::VertexID(iCell),
                {}, ///< edges of the adjacency graph are added automatically during edge insertion
                Mesh::Vertex::Data {
                    .diffusivity = 1.0,
                    .axes = axes,
                    .spatialTransform = maths::ScaleTranslateTransform<Scalar,Dimension>(transformed.begin(), transformed.end()),
                    .pAnsatzSpace = pAnsatzSpace,
                    .pAnsatzDerivatives = pAnsatzDerivatives
                }
            ));

            // Insert the current cell's connections to other cells already in the mesh.
            // The rule here is that cells with lower manhattan distance from the origin
            // are sources, while those with a higher norm are targets.
            if (iCellRow) {
                const Size iSourceCell = iCell - (nodesPerDirection - 1);
                const Size iTargetCell = iCell;
                BoundaryID source = iCellRow % 2 ? BoundaryID("+x") : BoundaryID("-x");
                BoundaryID target = iCellRow % 2 ? BoundaryID("+x") : BoundaryID("-x");
                mesh.insert(Mesh::Edge(
                    Mesh::EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    std::make_pair(source, target)
                ));
            } // if iCellRow

            if (iCellColumn) {
                const Size iSourceCell = iCell - 1ul;
                const Size iTargetCell = iCell;
                BoundaryID source = iCellColumn % 2 ? BoundaryID("+y") : BoundaryID("-y");
                BoundaryID target = iCellColumn % 2 ? BoundaryID("+y") : BoundaryID("-y");
                mesh.insert(Mesh::Edge(
                    Mesh::EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    std::make_pair(source, target)
                ));
            } // if iCellColumn
        } // for iCellColumn in range(nodesPerDirection -1)
    } // for iCellRow in range(nodesPerDirection - 1)

    Assembler assembler;
    assembler.addGraph(mesh,
                       [pAnsatzSpace]([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) -> std::size_t {
                            return pAnsatzSpace->size();
                        },
                       [&ansatzMap, &mesh](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {
                            const auto sourceAxes = mesh.find(rEdge.source()).value().data().axes;
                            const auto targetAxes = mesh.find(rEdge.target()).value().data().axes;
                            ansatzMap.getPairs(OrientedBoundary<Dimension>(sourceAxes, rEdge.data().first),
                                               OrientedBoundary<Dimension>(targetAxes, rEdge.data().second),
                                               it);
                        });

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<double> nonzeros;
    assembler.makeCSRMatrix(rowCount, columnCount, rowExtents, columnIndices, nonzeros);
    DynamicArray<Scalar> rhs(rowCount, 0.0);

    // Compute element contributions and assemble them into the matrix
    {
        const Quadrature<Scalar,Dimension> quadrature(GaussLegendreQuadrature<Scalar>(/*integrationOrder=*/2));
        DynamicArray<Scalar> derivativeBuffer(pAnsatzDerivatives->size());
        DynamicArray<Scalar> integrandBuffer(pAnsatzSpace->size() * pAnsatzSpace->size());
        DynamicArray<Scalar> productBuffer(integrandBuffer.size());

        for (Ref<const Mesh::Vertex> rCell : mesh.vertices()) {
            const auto jacobian = rCell.data().spatialTransform.makeDerivative();

            const auto localIntegrand = maths::makeTransformedIntegrand(
                maths::LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(rCell.data().diffusivity,
                                                                             rCell.data().pAnsatzDerivatives,
                                                                             {derivativeBuffer.data(), derivativeBuffer.size()}),
                jacobian
            );
            quadrature.evaluate(localIntegrand, integrandBuffer.data());

            {
                const auto keys = assembler.keys();
                CIE_TEST_REQUIRE(std::find(keys.begin(), keys.end(), rCell.id()) != keys.end());
            }
            const auto& rGlobalDofIndices = assembler[rCell.id()];
            const unsigned localSystemSize = rCell.data().pAnsatzSpace->size();

            for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
                for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
                    CIE_TEST_REQUIRE(iLocalRow < rGlobalDofIndices.size());
                    CIE_TEST_REQUIRE(iLocalColumn < rGlobalDofIndices.size());

                    const auto iRowBegin = rowExtents[rGlobalDofIndices[iLocalRow]];
                    const auto iRowEnd = rowExtents[rGlobalDofIndices[iLocalRow] + 1];
                    const auto itColumnIndex = std::find(columnIndices.begin() + iRowBegin,
                                                         columnIndices.begin() + iRowEnd,
                                                         rGlobalDofIndices[iLocalColumn]);
                    CIE_OUT_OF_RANGE_CHECK(itColumnIndex != columnIndices.begin() + iRowEnd)
                    const auto iEntry = std::distance(columnIndices.begin(),
                                                      itColumnIndex);
                    nonzeros[iEntry] += integrandBuffer[iLocalRow * localSystemSize + iLocalColumn];
                } // for iLocalColumn in range(ansatzBuffer.size)
            } // for iLocalRow in range(ansatzBuffer.size)
        } // for rCell in mesh.vertices
    }

    // Find DoFs to constrain:
    // 0) u(0, 0) = 0
    // 1) u(1, 0) = 1
    // 2) u(0, 1) = 2
    // 3) u(1, 1) = 3
    {
        StaticArray<std::optional<std::size_t>,4> iConstrainedDofs;
        utils::Comparison<Scalar> comparison(1e-8, 1e-6);

        for (const auto& rCell : mesh.vertices()) {
            // The DoFs we're looking for are nodal DoFs, meaning we need to find
            // a cell that has a vertex at a corner, and find the DoF within whose
            // basis function evaluates to 1 at that location.
            std::optional<std::size_t> maybeLocalCornerIndex;
            std::optional<std::size_t> maybeGlobalCornerIndex;

            for (const unsigned iLocalCorner : std::ranges::views::iota(0u, intPow(2u, Dimension))) {
                StaticArray<Scalar,Dimension> localCoordinates, globalCoordinates;
                localCoordinates[0] = (iLocalCorner & 1u) ? 1.0 : -1.0;
                localCoordinates[1] = (iLocalCorner & 2u) ? 1.0 : -1.0;
                rCell.data().spatialTransform.evaluate(localCoordinates.data(),
                                                       localCoordinates.data() + localCoordinates.size(),
                                                       globalCoordinates.data());

                // Loop over global corners and see whether any coincide with the transformed point.
                for (const unsigned iGlobalCorner : std::ranges::views::iota(0u, intPow(2u, Dimension))) {
                    StaticArray<Scalar,Dimension> referenceCoordinates;
                    referenceCoordinates[0] = Scalar(iGlobalCorner & 1u);
                    referenceCoordinates[1] = Scalar((iGlobalCorner & 2u) >> 1);
                    if (comparison.equal(globalCoordinates[0], referenceCoordinates[0]) && comparison.equal(globalCoordinates[1], referenceCoordinates[1])) {
                        maybeLocalCornerIndex = iLocalCorner;
                        maybeGlobalCornerIndex = iGlobalCorner;
                        break;
                    }
                } //for iGlobalCorner in range(4)

                if (maybeLocalCornerIndex.has_value()) {
                    break;
                }
            } // for iLocalCorner in range(4)

            // Found a cell that has a vertex at a corner.
            if (maybeLocalCornerIndex.has_value()) {
                StaticArray<Scalar,Dimension> localCoordinates;
                localCoordinates[0] = (maybeLocalCornerIndex.value() & 1u) ? 1.0 : -1.0;
                localCoordinates[1] = (maybeLocalCornerIndex.value() & 2u) ? 1.0 : -1.0;
                DynamicArray<Scalar> ansatzValues(rCell.data().pAnsatzSpace->size());
                rCell.data().pAnsatzSpace->evaluate(localCoordinates.data(),
                                                    localCoordinates.data() + localCoordinates.size(),
                                                    ansatzValues.data());
                for (unsigned iAnsatz=0u; iAnsatz<ansatzValues.size(); ++iAnsatz) {
                    if (comparison.equal(ansatzValues[iAnsatz], 1.0)) {
                        iConstrainedDofs[maybeGlobalCornerIndex.value()] = assembler[rCell.id()][iAnsatz];
                        break;
                    }
                }
            }
        } // for rCell in mesh.vertices

        // Barbaric imposition of dirichlet conditions.
        for (const auto maybeDofIndex : iConstrainedDofs) {
            CIE_TEST_REQUIRE(maybeDofIndex.has_value());
        }

        DynamicArray<std::pair<std::size_t,Scalar>> dirichletConditions;
        for (unsigned iDof=0u; iDof<iConstrainedDofs.size(); ++iDof) {
            dirichletConditions.emplace_back(iConstrainedDofs[iDof].value(), Scalar(iDof));
        }

        for (const auto& [iDof, value] : dirichletConditions) {
            for (std::size_t iRow=0ul; iRow<static_cast<std::size_t>(rowCount); ++iRow) {
                const std::size_t iEntryBegin = rowExtents[iRow];
                const std::size_t iEntryEnd   = rowExtents[iRow + 1];
                for (std::size_t iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                    const std::size_t iColumn = columnIndices[iEntry];
                    if (iRow == iDof) {
                        if (iRow == iColumn) nonzeros[iEntry] = 1.0;
                        else nonzeros[iEntry] = 0.0;
                    } else if (iColumn == iDof) {
                        rhs[iRow] -= value * nonzeros[iEntry];
                        nonzeros[iEntry] = 0.0;
                    }
                } // for iEntry in range(iEntryBegin, iEntryEnd)
            } // for iRow in range(rowCount)
            rhs[iDof] = value;
        }
    }

    // Solve the linear system.
    DynamicArray<Scalar> solution(rhs.size());
    {
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            rowCount,
            columnCount,
            nonzeros.size(),
            rowExtents.data(),
            columnIndices.data(),
            nonzeros.data());
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);

        Eigen::ConjugateGradient<EigenSparseMatrix,Eigen::Lower|Eigen::Upper,Eigen::DiagonalPreconditioner<Scalar>> solver;
        solver.setMaxIterations(int(1e3));
        solver.setTolerance(1e-6);

        solver.compute(lhsAdaptor);
        const auto x = solver.solve(rhsAdaptor);

        std::copy(x.begin(), x.end(), solution.begin());
    }

//    {
//        std::ofstream file("lhs.mm");
//        utils::io::MatrixMarket::Output io(file);
//        io(rowCount,
//           columnCount,
//           nonzeros.size(),
//           rowExtents.data(),
//           columnIndices.data(),
//           nonzeros.data());
//    }
//
//    {
//        std::ofstream file("rhs.mm");
//        utils::io::MatrixMarket::Output io(file);
//        io(rhs.data(), rhs.size());
//    }
//
//    {
//        std::ofstream file("u.mm");
//        utils::io::MatrixMarket::Output io(file);
//        io(solution.data(), solution.size());
//    }
//
//    {
//        std::ofstream file("mesh.gv");
//        io::Graphviz::Output io(file, io::Graphviz::Settings {.directed = false});
//        io(
//            mesh,
//            [&assembler](const Mesh::Vertex& rVertex) -> std::string {
//                std::stringstream label;
//                label << rVertex.id() << ":";
//                for (const auto& rItem : assembler[rVertex.id()]) {
//                    label << rItem << ',';
//                }
//                return label.str();
//            },
//            [](const Mesh::Edge& rEdge) -> std::string {
//                std::stringstream label;
//                label << rEdge.data().first << ":" << rEdge.data().second;
//                return label.str();
//            }
//        );
//    }

    // Postprocessing
    constexpr unsigned cellSamplesPerDirection = 2;
    DynamicArray<std::pair<StaticArray<Scalar,Dimension>,Scalar>> solutionSamples;
    solutionSamples.reserve(cellSamplesPerDirection * intPow(nodesPerDirection - 1u, Dimension));

    {
        DynamicArray<Scalar> ansatzBuffer(pAnsatzSpace->size());
        DynamicArray<StaticArray<Scalar,Dimension>> localSamplePoints;

        {
            StaticArray<unsigned,Dimension> samplePointState;
            std::fill(samplePointState.begin(), samplePointState.end(), 0u);
            do {
                StaticArray<Scalar,Dimension> localSamplePoint;
                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                    localSamplePoint[iDimension] = -1.0 + samplePointState[iDimension] * 2.0 / (cellSamplesPerDirection - 1.0);
                }
                localSamplePoints.emplace_back(localSamplePoint);
            } while (maths::OuterProduct<Dimension>::next(cellSamplesPerDirection, samplePointState.data()));
        }

        for (const auto& rCell : mesh.vertices()) {
            const auto& rGlobalIndices = assembler[rCell.id()];

            for (const auto& localCoordinates : localSamplePoints) {
                StaticArray<Scalar,Dimension> globalSamplePoint;
                rCell.data().spatialTransform.evaluate(localCoordinates.data(),
                                                       localCoordinates.data() + localCoordinates.size(),
                                                       globalSamplePoint.data());
                solutionSamples.emplace_back(globalSamplePoint, 0.0);
                ansatzBuffer.resize(rCell.data().pAnsatzSpace->size());
                rCell.data().pAnsatzSpace->evaluate(localCoordinates.data(),
                                                    localCoordinates.data() + localCoordinates.size(),
                                                    ansatzBuffer.data());

                for (unsigned iFunction=0u; iFunction<ansatzBuffer.size(); ++iFunction) {
                    solutionSamples.back().second += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }
            }
        }

        std::ofstream file("test_2d_solution.csv");
        for (const auto& [rSamplePoint, rValue] : solutionSamples) {
            for (const auto coordinate : rSamplePoint) file << coordinate << ',';
            file << rValue << '\n';
        }
    }
}


} // namespace cie::fem
