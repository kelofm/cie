// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"
#include "poisson2D/MeshData.hpp"
#include "poisson2D/CellData.hpp"
#include "poisson2D/BoundaryData.hpp"
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"
#include "packages/trees/inc//ContiguousSpaceTree.hpp"

// --- STL Includes ---
#include <packages/partitioning/inc/BoxBoundable.hpp>
#include <packages/primitives/inc/Object.hpp>
#include <ranges> // ranges::iota
#include <numbers> // std::numbers::pi


namespace cie::fem {


static_assert(::cie::concepts::SamplableGeometry<CellData>);
using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


/// @brief Data structure unique to the triangulated, immersed boundary cells.
using BoundaryCellData = maths::AffineEmbedding<Scalar,1u,Dimension>;


/// @brief Data structure unqiue to the triangulated, immersed boundary cell corners.
using BoundaryCornerData = Scalar;


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    BoundaryCornerData
>;


BVH makeBoundingVolumeHierarchy(Ref<Mesh> rMesh)
{
    auto logBlock = utils::LoggerSingleton::get().newBlock("make BVH");

    constexpr int targetLeafWidth = 5;
    constexpr int maxTreeDepth = 5;
    constexpr Scalar epsilon = 1e-3;

    geo::AABBoxNode<CellData> root;

    geo::AABBoxNode<CellData>::Point rootBase, rootLengths;
    std::fill(rootBase.begin(), rootBase.end(), -epsilon);
    std::fill(rootLengths.begin(), rootLengths.end(), 1.0 + 3.0 * epsilon);
    root = geo::AABBoxNode<CellData>(rootBase, rootLengths, nullptr);

    for (auto& rCell : rMesh.vertices()) {
        root.insert(&rCell.data());
    }

    root.partition(targetLeafWidth, maxTreeDepth);
    root.shrink();

    return BVH::flatten(
        root,
        [] (Ref<const CellData> rCellData) -> unsigned {return rCellData.id;},
        std::allocator<std::byte>());
}


BoundaryMesh generateBoundaryMesh(const unsigned resolution)
{
    using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;

    if (resolution < 3) CIE_THROW(Exception, "boundary resolution must be 3 or greater")

    BoundaryMesh boundary;

    DynamicArray<Point> corners;
    for (unsigned iSegment=0u; iSegment<resolution + 1; ++iSegment) {
        const Scalar arcParameter = iSegment * 2 * std::numbers::pi / resolution;
        corners.push_back(Point {
            boundaryRadius * std::cos(arcParameter) + 0.5,
            boundaryRadius * std::sin(arcParameter) + 0.5
        });
    } // for iSegment in range(resolution)

    for (unsigned iCorner=0u; iCorner<corners.size(); ++iCorner) {
        const StaticArray<Point,2> transformed {
            corners[iCorner],
            corners[(iCorner + 1) % corners.size()]
        };

        boundary.insert(BoundaryMesh::Vertex(
            boundary.vertices().size(),
            {},
            maths::AffineEmbedding<Scalar,1u,Dimension>(transformed)
        ));
    } // for iCorner in range(1, corners.size())

    for (unsigned iCorner=0; iCorner<corners.size(); ++iCorner) {
        boundary.insert(BoundaryMesh::Edge(
            iCorner,
            {iCorner, (iCorner + 1) % boundary.vertices().size()},
            static_cast<Scalar>(iCorner)
        ));
    } // for iCorner in range(corners.size())

    return boundary;
}


Ptr<const CellData> findContainingCell(Ref<const geo::AABBoxNode<CellData>> rBVH,
                                       Ref<const geo::AABBoxNode<CellData>::Point> rPoint)
{
    const auto pBVHNode = rBVH.find(rPoint);
    if (!pBVHNode) return nullptr;

    for (const auto pCellData : pBVHNode->contained()) {
        if (pCellData->at(rPoint)) {
            return pCellData;
        }
    }

    for (const auto pCellData : pBVHNode->intersected()) {
        if (pCellData->at(rPoint)) {
            return pCellData;
        }
    }
    return nullptr;
}


template <class TDofMap>
void addLHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        std::span<const int> rowExtents,
                        std::span<const int> columnIndices,
                        std::span<Scalar> entries)
{
    const unsigned localSystemSize = dofMap.size();
    for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
        for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
            const auto iRowBegin = rowExtents[dofMap[iLocalRow]];
            const auto iRowEnd = rowExtents[dofMap[iLocalRow] + 1];
            const auto itColumnIndex = std::lower_bound(columnIndices.begin() + iRowBegin,
                                                        columnIndices.begin() + iRowEnd,
                                                        dofMap[iLocalColumn]);
            CIE_OUT_OF_RANGE_CHECK(itColumnIndex != columnIndices.begin() + iRowEnd
                                && *itColumnIndex == dofMap[iLocalColumn]);
            const auto iEntry = std::distance(columnIndices.begin(),
                                            itColumnIndex);
            entries[iEntry] += contribution[iLocalRow * localSystemSize + iLocalColumn];
        } // for iLocalColumn in range(ansatzBuffer.size)
    } // for iLocalRow in range(ansatzBuffer.size)
}


template <class TDofMap>
void addRHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        std::span<Scalar> rhs)
{
    for (unsigned iComponent=0u; iComponent<contribution.size(); ++iComponent) {
        const auto iRow = dofMap[iComponent];
        rhs[iRow] += contribution[iComponent];
    }
}


struct DirichletBoundary : public maths::ExpressionTraits<Scalar>
{
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;

    void evaluate([[maybe_unused]] ConstSpan position, Span state) const noexcept {
        state[0] = position[0] + position[1];
    }

    unsigned size() const noexcept {
        return 1;
    }
}; // class DirichletBoundary


[[nodiscard]] DynamicArray<StaticArray<Scalar,2*Dimension+1>>
imposeBoundaryConditions(Ref<Mesh> rMesh,
                         Ref<const Assembler> rAssembler,
                         BVH::View bvh,
                         std::span<const CellData> contiguousCellData,
                         std::span<const int> rowExtents,
                         std::span<const int> columnIndices,
                         std::span<Scalar> entries,
                         std::span<Scalar> rhs,
                         Ref<const utils::ArgParse::Results> rArguments)
{
    auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

    DynamicArray<StaticArray<Scalar,2*Dimension+1>> boundarySegments; // {p0x, p0y, p1x, p1y, level}
    const unsigned boundaryResolution = rArguments.get<std::size_t>("boundary-resolution");
    const unsigned boundaryIntegrationOrder = rArguments.get<std::size_t>("boundary-integration-order");

    // Load the boundary mesh.
    const auto boundary = generateBoundaryMesh(boundaryResolution);

    using TreePrimitive = geo::Cube<1,Scalar>;
    using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

    const Quadrature<Scalar,1> lineQuadrature((GaussLegendreQuadrature<Scalar>(boundaryIntegrationOrder)));
    DynamicArray<Scalar> quadratureBuffer(  std::pow(rMesh.data().ansatzSpaces.front().size(),Dimension)
                                          + rMesh.data().ansatzSpaces.front().size());
    DynamicArray<Scalar> integrandBuffer(  rMesh.data().ansatzSpaces.front().size()
                                         + Dimension
                                         + 1);

    DirichletBoundary dirichletBoundary;
    Scalar boundaryLength = 0.0;

    for (const auto& rBoundaryCell : boundary.vertices()) {
        Tree tree(Tree::Point {-1.0}, 2.0);

        // Define a functor detecting cell boundaries.
        const auto boundaryVisitor = [bvh, &tree, &rBoundaryCell,
                                      &lineQuadrature, &integrandBuffer, &quadratureBuffer, &rMesh,
                                      &rowExtents, &columnIndices, &entries, &rAssembler,
                                      &rhs, &dirichletBoundary, &boundaryLength,
                                      &boundarySegments, &contiguousCellData, &rArguments] (
            Ref<const Tree::Node> rNode,
            unsigned level) -> bool {

            if (level < minBoundaryTreeDepth) return true;
            if (maxBoundaryTreeDepth < level) return false;

            Scalar base, edgeLength;
            tree.getNodeGeometry(rNode, &base, &edgeLength);

            // Define sample points in the boundary cell's local space.
            const StaticArray<Scalar,1>
                localBase {base},
                localOpposite {base + edgeLength};

            // Transform sample points to the global coordinate space.
            geo::Traits<Dimension,Scalar>::Point globalBase, globalOpposite;
            rBoundaryCell.data().evaluate(localBase, globalBase);
            rBoundaryCell.data().evaluate(localOpposite, globalOpposite);

            // Check whether the two endpoints are in different cells.
            const auto iMaybeBaseCell = bvh.find(
                std::span<const Scalar,Dimension>(globalBase.data(), Dimension),
                contiguousCellData
            );
            const auto iMaybeOppositeCell = bvh.find(
                std::span<const Scalar,Dimension>(globalOpposite.data(), Dimension),
                contiguousCellData
            );

            // Integrate if both endpoints lie in the same cell.
            if (iMaybeBaseCell != contiguousCellData.size() && iMaybeOppositeCell != contiguousCellData.size() && iMaybeBaseCell == iMaybeOppositeCell) {
                const Scalar segmentNorm =   std::pow(globalOpposite[0] - globalBase[0], static_cast<Scalar>(2))
                                           + std::pow(globalOpposite[1] - globalBase[1], static_cast<Scalar>(2));

                if (minBoundarySegmentNorm < segmentNorm) {
                    Ref<const CellData> rCell = contiguousCellData[iMaybeBaseCell];
                    const auto& rAnsatzSpace = rMesh.data().ansatzSpaces[rCell.iAnsatz];

                    StaticArray<maths::AffineEmbedding<Scalar,1,Dimension>::OutPoint,2> globalCorners;
                    globalCorners[0][0] = globalBase[0];
                    globalCorners[0][1] = globalBase[1];
                    globalCorners[1][0] = globalOpposite[0];
                    globalCorners[1][1] = globalOpposite[1];
                    const maths::AffineEmbedding<Scalar,1,Dimension> segmentTransform(globalCorners);

                    const auto integrand = makeTransformedIntegrand(
                        makeDirichletPenaltyIntegrand(dirichletBoundary,
                                                      /*penalty=*/rArguments.get<double>("penalty-factor"),
                                                      rAnsatzSpace,
                                                      segmentTransform,
                                                      std::span<Scalar>(integrandBuffer)),
                        segmentTransform.makeDerivative());
                    lineQuadrature.evaluate(integrand, quadratureBuffer);
                    const auto& rGlobalDofIndices = rAssembler[rCell.id];

                    addLHSContribution({quadratureBuffer.data(),
                                           static_cast<std::size_t>(std::pow(rAnsatzSpace.size(),Dimension))},
                                       rGlobalDofIndices,
                                       rowExtents,
                                       columnIndices,
                                       entries);

                    addRHSContribution({quadratureBuffer.data() + static_cast<std::size_t>(std::pow(rAnsatzSpace.size(),Dimension)),
                                           quadratureBuffer.data() + quadratureBuffer.size()},
                                       rGlobalDofIndices,
                                       rhs);

                    // Log debug and output info.
                    boundaryLength += std::sqrt(segmentNorm);
                    decltype(boundarySegments)::value_type segment;
                    std::copy_n(globalCorners[0].data(), Dimension, segment.data());
                    std::copy_n(globalCorners[1].data(), Dimension, segment.data() + Dimension);
                    segment.back() = level;
                    boundarySegments.push_back(segment);
                }
            } // if both endpoints lie in the same cell

            return iMaybeBaseCell != iMaybeOppositeCell &&
                   (iMaybeBaseCell != contiguousCellData.size() || iMaybeOppositeCell != contiguousCellData.size());
        }; // boundaryVisitor

        // Construct a binary tree that detects intersections between
        // Cell boundaries and the current boundary cell.
        tree.scan(boundaryVisitor);
    } // for rBoundaryCell in boundary.vertices()

    return boundarySegments;
}


/// @brief 2D system test.
int main(Ref<const utils::ArgParse::Results> rArguments) {
    Mesh mesh;
    mp::ThreadPoolBase threads;

    const unsigned integrationOrder = rArguments.get<std::size_t>("i");
    const unsigned postprocessResolution = rArguments.get<std::size_t>("scatter-resolution");

    // Fill the mesh with cells and boundaries.
    generateMesh(mesh, rArguments);

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("write graphml");
        io::GraphML::GraphML::Output("test_2d.graphml")(mesh);
    }

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent cells, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const StaticArray<Scalar,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    const auto ansatzMap = makeAnsatzMap(mesh.data().ansatzSpaces.front(),
                                         samples,
                                         utils::Comparison<Scalar>(/*absoluteTolerance =*/ 1e-8,
                                                                   /*relativeTolerance =*/ 1e-6));

    // Build a factory detects the topology of the mesh and issues indices to DoFs.
    Assembler assembler;

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("parse mesh topology");
        assembler.addGraph(mesh,
                           [&mesh]([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) -> std::size_t {
                                const Ansatz& rAnsatz = mesh.data().ansatzSpaces[rVertex.data().iAnsatz];
                                return rAnsatz.size();
                           },
                           [&ansatzMap, &mesh](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {
                                const auto sourceAxes = mesh.find(rEdge.source()).value().data().axes;
                                const auto targetAxes = mesh.find(rEdge.target()).value().data().axes;
                                ansatzMap.getPairs(OrientedBoundary<Dimension>(sourceAxes, rEdge.data().boundary),
                                                   OrientedBoundary<Dimension>(targetAxes, rEdge.data().boundary),
                                                   it);
                           });
    } // parse mesh topology

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<double> entries;
    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("compute sparsity pattern");
        assembler.makeCSRMatrix(rowCount, columnCount, rowExtents, columnIndices, entries);
    }
    DynamicArray<Scalar> rhs(rowCount, 0.0);

    // Compute element contributions and assemble them into the matrix
    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");

        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(integrationOrder)));
        DynamicArray<Scalar> derivativeBuffer(mesh.data().ansatzDerivatives.front().size());
        DynamicArray<Scalar> integrandBuffer(std::pow(mesh.data().ansatzSpaces.front().size(), Dimension));
        DynamicArray<Scalar> productBuffer(integrandBuffer.size());

        for (Ref<const Mesh::Vertex> rCell : mesh.vertices()) {
            auto& rAnsatzDerivatives  = mesh.data().ansatzDerivatives[rCell.data().iAnsatz];
            const auto jacobian = rCell.data().spatialTransform.makeDerivative();

            const auto localIntegrand = makeTransformedIntegrand(
                LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(rCell.data().diffusivity,
                                                                      rAnsatzDerivatives,
                                                                      {derivativeBuffer.data(), derivativeBuffer.size()}),
                jacobian
            );
            quadrature.evaluate(localIntegrand, integrandBuffer);
            const auto& rGlobalDofIndices = assembler[rCell.id()];
            addLHSContribution(integrandBuffer,
                               rGlobalDofIndices,
                               rowExtents,
                               columnIndices,
                               entries);
        } // for rCell in mesh.vertices
    } // integrate

    // Construct a bounding volume hierarchy over the cells to accelerate
    // point membership tests. Running on accelerator devices also requires
    // - the cell data to be available in a contiguous array
    // - the cell data to be self contained (no pointers and heap storage)
    auto bvh = makeBoundingVolumeHierarchy(mesh);
    const auto bvhView = bvh.makeView();
    DynamicArray<CellData> contiguousCellData(mesh.vertices().size());
    std::ranges::transform(
        mesh.vertices(),
        contiguousCellData.data(),
        [](Ref<const Mesh::Vertex> rCell){return rCell.data();});

    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        assembler,
        bvhView,
        contiguousCellData,
        rowExtents,
        columnIndices,
        entries,
        rhs,
        rArguments);

    // Solve the linear system.
    DynamicArray<Scalar> solution(rhs.size());
    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("solve linear system");
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            rowCount,
            columnCount,
            entries.size(),
            rowExtents.data(),
            columnIndices.data(),
            entries.data());
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);

        Eigen::ConjugateGradient<
            EigenSparseMatrix,
            Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<Scalar>
        > solver;
        solver.setMaxIterations(int(5e3));
        solver.setTolerance(1e-6);

        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);

        std::cout << solver.iterations() << " iterations "
                  << solver.error()      << " residual\n";
    } // solve

    // XDMF output.
    {
        // Collect output for the bounding volume hierarchy.
        DynamicArray<StaticArray<Scalar,Dimension*intPow(2,Dimension)>> boundingVolumes; // {p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y}
        bvh.visit([&boundingVolumes](const auto& rBox) -> bool {
            constexpr unsigned cornerCount = intPow(2,Dimension);
            boundingVolumes.emplace_back();
            geo::Box<Dimension,Scalar>::makeCorners(
                rBox.base(),
                rBox.lengths(),
                std::span<Scalar,cornerCount*Dimension> (boundingVolumes.back().data(), cornerCount * Dimension)
            );
            return true;
        });

        // Collect state samples.
        struct Sample {
            StaticArray<CellData::GlobalCoordinate,Dimension>   position;
            Scalar                                              state;
            unsigned                                            cellID;
        }; // struct Sample
        DynamicArray<Sample> samples(postprocessResolution * postprocessResolution);

        {
            auto logBlock = utils::LoggerSingleton::get().newBlock("scatter postprocess");
            constexpr Scalar epsilon = 1e-10;
            const Scalar postprocessDelta  = (1.0 - 2 * epsilon) / (postprocessResolution - 1);
            DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.front().size());

            mp::ParallelFor<unsigned>(threads).firstPrivate(DynamicArray<Scalar>(), mesh.data().ansatzSpaces)(
                intPow(postprocessResolution, 2),
                [&samples, &solution, &assembler, bvhView, &contiguousCellData, postprocessResolution, postprocessDelta](
                        const unsigned iSample,
                        Ref<DynamicArray<Scalar>> rAnsatzBuffer,
                        Ref<const DynamicArray<Ansatz>> rAnsatzSpaces) -> void {
                    const unsigned iSampleY = iSample / postprocessResolution;
                    const unsigned iSampleX = iSample % postprocessResolution;
                    auto& rSample = samples[iSample];
                    rSample.position = {epsilon + iSampleX * postprocessDelta,
                                        epsilon + iSampleY * postprocessDelta};

                    // Find which cell the global point lies in.
                    const auto iMaybeCellData = bvhView.find(
                        Kernel<Dimension,Scalar>::decayView(rSample.position),
                        std::span<const CellData>(contiguousCellData)
                    );

                    if (iMaybeCellData != contiguousCellData.size()) {
                        Ref<const CellData> rCellData = contiguousCellData[iMaybeCellData];
                        rSample.cellID = rCellData.id;

                        // Compute sample point in the cell's local space.
                        StaticArray<CellData::LocalCoordinate,Dimension> localSamplePoint;
                        rCellData.transform(
                            Kernel<Dimension,Scalar>::view(rSample.position),
                            Kernel<Dimension,Scalar>::view(localSamplePoint));

                        // Evaluate the cell's ansatz functions at the local sample point.
                        const auto& rAnsatzSpace = rAnsatzSpaces[rCellData.iAnsatz];
                        rAnsatzBuffer.resize(rAnsatzSpace.size());
                        rAnsatzSpace.evaluate(Kernel<Dimension,Scalar>::decayView(localSamplePoint), rAnsatzBuffer);

                        // Find the entries of the cell's DoFs in the global state vector.
                        const auto& rGlobalIndices = assembler[rCellData.id];

                        // Compute state as an indirect inner product of the solution vector
                        // and the ansatz function values at the local corner coordinates.
                        rSample.state = 0;
                        for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                            rSample.state += solution[rGlobalIndices[iFunction]] * rAnsatzBuffer[iFunction];
                        } // for iFunction in range(rGlobalIndices.size())
                    } else {
                        // Could not find the cell that contains the current sample point.
                        rSample.state = NAN;
                    }
                });
        }

        // Write XDMF.
        std::ofstream xdmf("test_2d.xdmf");
        xdmf << R"(
<Xdmf Version="3.0">
    <Domain>
        <Grid name="root" GridType="Collection" CollectionType="Spatial">

            <Grid Name="boundary" GridType="Uniform">
                <Topology TopologyType="Polyline" NumberOfElements=")" << boundarySegments.size() << R"(" NodesPerElement="2">
                    <DataItem Format="XML" Dimensions=")" << boundarySegments.size() << R"( 2">
                        )"; {
                            for (unsigned iSegment=0u; iSegment<boundarySegments.size(); ++iSegment) {
                                xdmf << 2 * iSegment << " " << 2 * iSegment + 1 << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format="XML" Dimensions=")" << 2 * boundarySegments.size() << R"( 3">
                        )"; {
                            for (const auto& rSegment : boundarySegments) {
                                xdmf << rSegment[0] << " " << rSegment[1] << " 0" << "\n                        "
                                     << rSegment[2] << " " << rSegment[3] << " 0" << "\n                        ";
                            }
                        };
                        xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="level" Center="Cell" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << boundarySegments.size() << R"( 1">
                        )"; {
                            for (const auto& rSegment : boundarySegments) {
                                xdmf << rSegment.back() << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>

            <Grid Name="bvh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << boundingVolumes.size() << R"(" NodesPerElement="4">
                    <DataItem Format="XML" Dimensions=")" << boundingVolumes.size() << R"( 4">
                        )"; {
                            for (unsigned iBoundingVolume=0u; iBoundingVolume<boundingVolumes.size(); ++iBoundingVolume) {
                                const unsigned iPointBegin = 4 * iBoundingVolume;
                                xdmf << iPointBegin << " "
                                     << iPointBegin + 1 << " "
                                     << iPointBegin + 3 << " "
                                     << iPointBegin + 2
                                     << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format="XML" Dimensions=")" << boundingVolumes.size() * 4 << R"( 3">
                        )"; {
                            for (const auto& rBoundingVolume : boundingVolumes) {
                                for (unsigned iCorner=0u; iCorner<4u; ++iCorner) {
                                    for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension)
                                        xdmf << rBoundingVolume[Dimension * iCorner + iDimension] << " ";
                                    xdmf << 0 << "\n                        ";
                                }
                            }
                        } xdmf << R"(
                    </DataItem>
                </Geometry>
            </Grid>

            <Grid Name="mesh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << mesh.vertices().size() << R"(" NodesPerElement="4">
                    <DataItem Format="XML" Dimensions=")" << mesh.vertices().size() << R"( 4">
                        )"; {
                            for (unsigned iCell=0u; iCell<mesh.vertices().size(); ++iCell) {
                                const unsigned i = 4 * iCell;
                                xdmf << i << " " << i + 1 << " " << i + 3 << " " << i + 2 << "\n                        ";
                            } // for iCell in range(mesh.vertices().size())
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions=")" << 4 * mesh.vertices().size() << R"( 3">
                        )"; {
                            StaticArray<StaticArray<CellData::LocalCoordinate,Dimension>,intPow(2,Dimension)> localCorners {
                                {-1.0, -1.0},
                                { 1.0, -1.0},
                                {-1.0,  1.0},
                                { 1.0,  1.0}
                            };

                            for (const auto& rCell : mesh.vertices()) {
                                for (const auto& rLocalPoint : localCorners) {
                                    StaticArray<CellData::GlobalCoordinate,Dimension> globalPoint;
                                    rCell.data().transform(
                                        Kernel<Dimension,Scalar>::view(rLocalPoint),
                                        Kernel<Dimension,Scalar>::view(globalPoint));
                                    for (auto c : globalPoint) xdmf << c << " ";
                                    xdmf << "0\n                        ";
                                } // for rLocalPoint : localCorners
                            } // for rCell in mesh.vertices()
                        } xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << 4 * mesh.vertices().size() << R"( 1">
                        )"; {
                            StaticArray<StaticArray<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
                                {-1.0, -1.0},
                                { 1.0, -1.0},
                                {-1.0,  1.0},
                                { 1.0,  1.0}
                            };
                            DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.front().size());

                            for (const auto& rCell : mesh.vertices()) {
                                const auto& rGlobalIndices = assembler[rCell.id()];
                                const auto& rAnsatzSpace = mesh.data().ansatzSpaces[rCell.data().iAnsatz];
                                ansatzBuffer.resize(rAnsatzSpace.size());

                                for (const auto& rLocalPoint : localCorners) {
                                    // Compute ansatz function values at the current local corner.
                                    rAnsatzSpace.evaluate(rLocalPoint, ansatzBuffer);

                                    // Compute state as an indirect inner product of the solution vector
                                    // and the ansatz function values at the local corner coordinates.
                                    Scalar state = 0;
                                    for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                                        state += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                                    } // for iFunction in range(rGlobalIndices.size())

                                    xdmf << state << "\n                        ";
                                } // for rLocalPoint : localCorners
                            } // for rCell in mesh.vertices()
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>

            <Grid Name="sample" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << intPow(postprocessResolution - 1, 2) << R"(" NodePerElement="4">
                    <DataItem Format="XML" Dimensions=")" << intPow(postprocessResolution - 1, 2) << R"( 4">
                        )"; {
                            for (unsigned iSampleCellY=0u; iSampleCellY<postprocessResolution - 1; ++iSampleCellY) {
                                for (unsigned iSampleCellX=0u; iSampleCellX<postprocessResolution - 1; ++iSampleCellX) {
                                    const unsigned iBase = iSampleCellY * postprocessResolution + iSampleCellX;
                                    xdmf << iBase << " " << iBase + 1 << " " << iBase + postprocessResolution +1 << " " << iBase + postprocessResolution << "\n                        ";
                                } // for iSampleCellX in range(proceprocessResolution)
                            } // for iSampleCellY in range(proceprocessResolution)
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 3">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.position[0] << " "
                                     << rSample.position[1]
                                     << " 0\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 1">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.state << "\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Attribute>

                <Attribute Name="cellID" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 1">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.cellID << "\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>
        </Grid>
    </Domain>
</Xdmf>
)";
    }

    return 0;
}


} // namespace cie::fem

int main(int argc, const char** argv) {
    cie::utils::ArgParse parser("2D Poisson Example");
    parser
        .addKeyword(
            {"-r", "--resolution"},
            cie::utils::ArgParse::DefaultValue {"31"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Resolution (number of nodes per direction.)")
        .addKeyword(
            {"-p", "--polynomial-order"},
            cie::utils::ArgParse::DefaultValue {"5"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Polynomial order of the ansatz space.")
        .addKeyword(
            {"-i", "--integration-order"},
            cie::utils::ArgParse::DefaultValue {"6"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Integration order for volumetric integrands.")
        .addKeyword(
            {"--boundary-integration-order"},
            cie::utils::ArgParse::DefaultValue {"6"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Integration order for boundary integrands.")
        .addKeyword(
            {"--boundary-resolution"},
            cie::utils::ArgParse::DefaultValue {"20"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of nodes discretizing the boundary circle.")
        .addKeyword(
            {"--scatter-resolution"},
            cie::utils::ArgParse::DefaultValue {"100"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of sample points per direction for scatter post-processing.")
        .addKeyword(
            {"--penalty-factor"},
            cie::utils::ArgParse::DefaultValue {"1e-2"},
            cie::utils::ArgParse::validatorFactory<double>(),
            "Penalty value for the weak imposition of Dirichlet boundary conditions.")
        .addFlag(
            {"-h", "--help"},
            "Print this help and exit.")
        ;

    cie::utils::ArgParse::Results arguments;
    try {
        arguments = parser.parseArguments(argc - 1, argv + 1);
    } catch (cie::Exception rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }

    if (arguments.get<bool>("h")) {
        parser.help(std::cout);
        return 0;
    }

    try {
        return cie::fem::main(arguments);
    } catch (cie::Exception& rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }
}
