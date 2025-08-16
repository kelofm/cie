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
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"

// --- STL Includes ---
#include <ranges> // ranges::iota


namespace cie::fem {


using Scalar = double;


template <class TAnsatzSpace>
struct CellData
{
    using SpatialTransform = maths::ScaleTranslateTransform<Scalar,1>;

    Scalar                                              diffusivity;
    OrientedAxes<TAnsatzSpace::Dimension>               axes;
    SpatialTransform                                    spatialTransform;
    std::shared_ptr<TAnsatzSpace>                       pAnsatzSpace;
    std::shared_ptr<typename TAnsatzSpace::Derivative>  pAnsatzDerivatives;
}; // struct CellData


/** @brief 1D system test.
 *  @details This test models linear steady-state heat convection on a 1D bar
 *           with Dirichlet conditions on both boundaries.
 *
 *           The mesh conforms to the boundaries and consists of 10 1D elements
 *           with linear basis functions.
 *
 *           Mesh:
 *           @code
 *                    +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+
 *           [u(0)=1] | 0 |-0-| 1 |-1-| 2 |-2-| 3 |-3-| 4 |-4-| 5 |-5-| 6 |-6-| 7 |-7-| 8 |-8-| 9 | [u(L)=0]
 *                    +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+
 *           @endcode
 */
CIE_TEST_CASE("1D", "[systemTests]")
{
    CIE_TEST_CASE_INIT("1D")
    constexpr unsigned Dimension = 1;
    constexpr Size nodeCount = 10;

    // Construct a linear 1D ansatz space
    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;
    auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
        //Basis({ 0.0, -0.5,  0.5}),
        //Basis({ 0.0,  0.5,  0.5}),
        //Basis({ 1.0,  0.0, -1.0})
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
    const Scalar elementSize = 1.0 / (nodeCount - 1);
    for (Size iCell : std::ranges::views::iota(0ul, nodeCount - 1)) {
        StaticArray<StaticArray<Scalar,Dimension>,2> transformed;
        OrientedAxes<Dimension> axes;

        if (iCell % 2) {
            axes[0] = "-x";
            transformed[0].front() = (iCell + 1.0) * elementSize;
            transformed[1].front() = iCell * elementSize;
        } else {
            axes[0] = "+x";
            transformed[0].front() = iCell * elementSize;
            transformed[1].front() = (iCell + 1.0) * elementSize;
        }

        // Insert the cell into the adjacency graph (mesh) as a vertex
        mesh.insert(Mesh::Vertex(
            VertexID(iCell),
            {}, ///< edges of the adjacency graph are added automatically during edge insertion
            Mesh::Vertex::Data {
                .diffusivity = 1.0,
                .axes = axes,
                .spatialTransform = maths::ScaleTranslateTransform<Scalar,Dimension>(transformed.begin(), transformed.end()),
                .pAnsatzSpace = pAnsatzSpace,
                .pAnsatzDerivatives = pAnsatzDerivatives
            }
        ));
    }

    // Insert shared element boundaries into the adjacency graph as edges
    for (Size iBoundary : std::ranges::views::iota(0ul, nodeCount - 2)) {
        const Size iLeftCell = iBoundary;
        const Size iRightCell = iBoundary + 1;

        BoundaryID left("+x"), right("+x");
        if (iBoundary % 2) {
            left = "-x";
            right = "-x";
        }

        mesh.insert(Mesh::Edge(
            iBoundary,
            {iLeftCell, iRightCell},
            std::make_pair(left, right) // <== all bundaries connect cells on the left to cells on the right
        ));
    }

    Assembler assembler;
    assembler.addGraph(mesh,
                       [pAnsatzSpace]([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) -> std::size_t {
                            return pAnsatzSpace->size();
                        },
                       [&ansatzMap, &mesh](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {
                            const auto sourceAxes = mesh.find(rEdge.source()).value().data().axes;
                            const auto targetAxes = mesh.find(rEdge.target()).value().data().axes;
                            ansatzMap.getPairs(OrientedBoundary<1>(sourceAxes, rEdge.data().first),
                                               OrientedBoundary<1>(targetAxes, rEdge.data().second),
                                               it);
                        });

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<double> nonzeros;
    assembler.makeCSRMatrix(rowCount, columnCount, rowExtents, columnIndices, nonzeros);

    // Compute element contributions and assemble them into the matrix
    {
        const Quadrature<Scalar,1> quadrature(GaussLegendreQuadrature<Scalar>(/*integrationOrder=*/2));
        DynamicArray<Scalar> derivativeBuffer(pAnsatzDerivatives->size());
        DynamicArray<Scalar> integrandBuffer(pAnsatzSpace->size() * pAnsatzSpace->size());
        DynamicArray<Scalar> productBuffer(integrandBuffer.size());

        for (Ref<const Mesh::Vertex> rCell : mesh.vertices()) {
            const auto jacobian = rCell.data().spatialTransform.makeDerivative();

            const auto localIntegrand = maths::makeLambdaExpression<Scalar>(
                [&jacobian, &derivativeBuffer, &rCell, &productBuffer]
                        (std::span<const Scalar> in, std::span<Scalar> out) -> void {
                    const Scalar jacobianDeterminant = jacobian.evaluateDeterminant(in);
                    rCell.data().pAnsatzDerivatives->evaluate(in, derivativeBuffer);

                    using MatrixAdaptor = Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>;
                    MatrixAdaptor derivativeAdaptor(derivativeBuffer.data(),
                                                    Dimension,
                                                    rCell.data().pAnsatzSpace->size());
                    MatrixAdaptor productAdaptor(productBuffer.data(),
                                                 rCell.data().pAnsatzSpace->size(),
                                                 rCell.data().pAnsatzSpace->size());
                    productAdaptor = derivativeAdaptor.transpose() * derivativeAdaptor;

                    auto itOut = out.data();
                    for (unsigned iLocalRow=0u; iLocalRow<productAdaptor.rows(); ++iLocalRow) {
                        for (unsigned iLocalColumn=0u; iLocalColumn<productAdaptor.cols(); ++iLocalColumn) {
                            *itOut++ += rCell.data().diffusivity * productAdaptor(iLocalRow, iLocalColumn) * std::abs(jacobianDeterminant);
                        } // for iLocalColumn in range(ansatzBuffer.size)
                    } // for iLocalRow in range(ansatzBuffer.size)
                },
                integrandBuffer.size());

            quadrature.evaluate(localIntegrand, integrandBuffer);

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

    // Find boundary DoFs
    std::optional<std::size_t> iLeftmostDof, iRightmostDof;
    {
        auto rLeftmostCell = mesh.findVertex(0);
        CIE_TEST_REQUIRE(rLeftmostCell.has_value());

        utils::Comparison<Scalar> comparison(1e-8, 1-6);
        DynamicArray<Scalar> ansatzBuffer(rLeftmostCell.value().data().pAnsatzSpace->size());
        StaticArray<Scalar,1> localCoordinates, globalCoordinates;

        localCoordinates.front() = -1.0;
        rLeftmostCell.value().data().spatialTransform.evaluate(localCoordinates, globalCoordinates);
        if (!comparison.equal(globalCoordinates.front(), 0.0)) {
            localCoordinates.front() = 1.0;
            rLeftmostCell.value().data().spatialTransform.evaluate(localCoordinates, globalCoordinates);
            CIE_TEST_REQUIRE(comparison.equal(globalCoordinates.front(), 0.0));
        }

        rLeftmostCell.value().data().pAnsatzSpace->evaluate(localCoordinates, ansatzBuffer);
        for (unsigned iAnsatz=0u; iAnsatz<ansatzBuffer.size(); ++iAnsatz) {
            if (comparison.equal(ansatzBuffer[iAnsatz], 1.0)) {
                iLeftmostDof = assembler[0][iAnsatz];
                break;
            }
        }
    }

    {
        auto rRightmostCell = mesh.findVertex(nodeCount - 2);
        CIE_TEST_REQUIRE(rRightmostCell.has_value());

        utils::Comparison<Scalar> comparison(1e-8, 1-6);
        DynamicArray<Scalar> ansatzBuffer(rRightmostCell.value().data().pAnsatzSpace->size());
        StaticArray<Scalar,1> localCoordinates, globalCoordinates;

        localCoordinates.front() = -1.0;
        rRightmostCell.value().data().spatialTransform.evaluate(localCoordinates, globalCoordinates);
        if (!comparison.equal(globalCoordinates.front(), 1.0)) {
            localCoordinates.front() = 1.0;
            rRightmostCell.value().data().spatialTransform.evaluate(localCoordinates, globalCoordinates);
            CIE_TEST_REQUIRE(comparison.equal(globalCoordinates.front(), 1.0));
        }

        rRightmostCell.value().data().pAnsatzSpace->evaluate(localCoordinates, ansatzBuffer);
        for (unsigned iAnsatz=0u; iAnsatz<ansatzBuffer.size(); ++iAnsatz) {
            if (comparison.equal(ansatzBuffer[iAnsatz], 1.0)) {
                iRightmostDof = assembler[nodeCount - 2][iAnsatz];
                break;
            }
        }
    }

    CIE_TEST_CHECK(iLeftmostDof.has_value());
    CIE_TEST_CHECK(iRightmostDof.has_value());

    // Barbaric dirichlet conditions.
    DynamicArray<Scalar> rhs(rowCount, 0.0);
    DynamicArray<std::pair<std::size_t,Scalar>> dirichletConditions;
    dirichletConditions.emplace_back(iLeftmostDof.value(), 1.0);
    dirichletConditions.emplace_back(iRightmostDof.value(), 0.0);

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

        Eigen::ConjugateGradient<EigenSparseMatrix,Eigen::Lower|Eigen::Upper,Eigen::IncompleteLUT<Scalar>> solver;
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

    // Postprocessing
    constexpr unsigned samplesPerCell = 5;
    DynamicArray<std::pair<StaticArray<Scalar,1>,Scalar>> solutionSamples;
    solutionSamples.reserve(samplesPerCell * (nodeCount - 1));

    {
        DynamicArray<Scalar> ansatzBuffer(pAnsatzSpace->size());
        DynamicArray<StaticArray<Scalar,1>> sampleCoordinates;

        for (unsigned iCoordinate=0u; iCoordinate<samplesPerCell; ++iCoordinate) {
            sampleCoordinates.emplace_back();
            sampleCoordinates.back().front() = -1.0 + iCoordinate * 2.0 / (samplesPerCell - 1.0);
        }

        for (const auto& rCell : mesh.vertices()) {
            const auto& rGlobalIndices = assembler[rCell.id()];

            for (const auto& localCoordinates : sampleCoordinates) {
                StaticArray<Scalar,1> globalCoordinates;
                rCell.data().spatialTransform.evaluate(localCoordinates, globalCoordinates);
                solutionSamples.emplace_back(globalCoordinates, 0.0);
                ansatzBuffer.resize(rCell.data().pAnsatzSpace->size());
                rCell.data().pAnsatzSpace->evaluate(localCoordinates, ansatzBuffer);

                for (unsigned iFunction=0u; iFunction<ansatzBuffer.size(); ++iFunction) {
                    solutionSamples.back().second += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }
            }
        }

        std::ofstream file("test_1d_solution.csv");
        for (const auto& [rCoordinates, rValue] : solutionSamples) {
            file << rCoordinates.front() << ',' << rValue << '\n';
        }
    }
}


} // namespace cie::fem
