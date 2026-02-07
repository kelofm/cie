#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/numeric/inc/KDTreeQuadraturePointFactory.hpp"
#include "packages/numeric/inc/CompositeDomain.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"

// --- STL Includes ---
#include <numbers> // std::numbers::pi


namespace cie::fem {


/// @brief Data structure unique to the triangulated, immersed boundary cells.
class BoundaryCellData : public CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension> {
public:
    using Base = CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension>;

    using Base::Base;

    BoundaryCellData(RightRef<typename Base::SpatialTransform> rEmbedding,
                     Ref<const CellData> rCell) noexcept
        : Base(
            VertexID(),
            Base::AnsatzSpaceID(),
            OrientedAxes<1>(),
            std::move(rEmbedding)),
         _pCell(&rCell)
    {}

    Ref<const CellData> getEmbeddingCell() const noexcept {
        return *_pCell;
    }

private:
    Ptr<const CellData> _pCell;
}; // class BoundaryCellData


class BoundaryMeshData {
public:
    BoundaryMeshData() noexcept = default;

    BoundaryMeshData(std::span<const QuadraturePoint<1,Scalar>> quadraturePointSet,
                     Ref<const MeshData> rMeshData)
        : _quadraturePointSet(quadraturePointSet.begin(), quadraturePointSet.end()),
          _pMeshData(&rMeshData)
    {}

    CachedQuadraturePointFactory<1,Scalar> makeQuadratureRule() const {
        return CachedQuadraturePointFactory<1,Scalar>(
            std::span<const QuadraturePoint<1,Scalar>>(
                _quadraturePointSet.data(),
                _quadraturePointSet.size()));
    }

private:
    std::vector<QuadraturePoint<1,Scalar>> _quadraturePointSet;

    Ptr<const MeshData> _pMeshData;
}; // class BoundaryMeshData


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    void,
    BoundaryMeshData
>;


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


struct DirichletBoundary : public maths::ExpressionTraits<Scalar> {
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;

    void evaluate([[maybe_unused]] ConstSpan position, Span state) const noexcept {
        state[0] = position[0] + position[1];
        //state[0] = 1.0;
    }

    unsigned size() const noexcept {
        return 1;
    }
}; // class DirichletBoundary


BVH makeBoundingVolumeHierarchy(Ref<Mesh> rMesh) {
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
        [] (Ref<const CellData> rCellData) -> unsigned {return rCellData.id();},
        std::allocator<std::byte>());
}


/// @brief Definition of a boundary segment in parametric space.
struct ParametricBoundarySegment {
    std::size_t iCell;
    Scalar segmentBegin;
    Scalar segmentEnd;
}; // struct ParametricBoundarySegment


void partitionBoundaryCell(Ref<maths::AffineEmbedding<Scalar,1u,Dimension>> rTransform,
                           BVH::View bvh,
                           std::span<const CellData> contiguousCellData,
                           unsigned minBoundaryTreeDepth,
                           unsigned maxBoundaryTreeDepth,
                           Ref<DynamicArray<ParametricBoundarySegment>> rOutput) {
    rOutput.clear();

    CIE_BEGIN_EXCEPTION_TRACING
    using TreePrimitive = geo::Cube<1,Scalar>;
    using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;
    Tree tree(Tree::Point {-1.0}, 2.0);

    const auto visitor = [&] (Ref<const Tree::Node> rNode, unsigned level) -> bool {
        // Skip processing if the minimum depth has not yet been reached.
        if (level < minBoundaryTreeDepth) return true;

        // Stop processing if the maximum depth is reached.
        if (maxBoundaryTreeDepth < level) return false;

        // Recover the node's geometry (line segment in parametric space).
        Scalar base, edgeLength;
        tree.getNodeGeometry(rNode, &base, &edgeLength);

        // Define sample points in the boundary cell's parametric space.
        const std::array<Scalar,1>
            parametricBase {base},
            parametricOpposite {base + edgeLength};

        // Transform sample points to the global coordinate space.
        std::array<Scalar,Dimension> physicalBase, physicalOpposite;
        rTransform.evaluate(parametricBase, physicalBase);
        rTransform.evaluate(parametricOpposite, physicalOpposite);

        // Check whether the two endpoints are in different cells.
        const std::size_t iBaseCell = bvh.find(
            std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                physicalBase.data()),
                Dimension),
            contiguousCellData);
        const std::size_t iOppositeCell = bvh.find(
            std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                physicalOpposite.data()),
                Dimension),
            contiguousCellData);

        // Extend the output if both endpoints lie in the same cell.
        if (iBaseCell != contiguousCellData.size() && iOppositeCell != contiguousCellData.size() && iBaseCell == iOppositeCell) {
            auto itSegment = std::lower_bound(
                rOutput.begin(),
                rOutput.end(),
                iBaseCell,
                [] (Ref<const ParametricBoundarySegment> rSegment, std::size_t iBaseCell) {
                    return iBaseCell < rSegment.iCell;
                });
            if (itSegment == rOutput.end()) {
                rOutput.emplace_back(ParametricBoundarySegment {
                    .iCell          = iBaseCell,
                    .segmentBegin   = base,
                    .segmentEnd     = base + edgeLength});
            } else if (itSegment->iCell != iBaseCell) {
                rOutput.insert(
                    itSegment,
                    ParametricBoundarySegment {
                        .iCell          = iBaseCell,
                        .segmentBegin   = base,
                        .segmentEnd     = base + edgeLength});
            } else {
                itSegment->segmentBegin = std::min(itSegment->segmentBegin, base             );
                itSegment->segmentEnd   = std::max(itSegment->segmentEnd,   base + edgeLength);
            }
        }

        return iBaseCell != iOppositeCell
            && (iBaseCell != contiguousCellData.size() || iOppositeCell != contiguousCellData.size());
    }; // visitor

    tree.scan(visitor);
    CIE_END_EXCEPTION_TRACING
}


BoundaryMesh generateBoundaryMesh(Ref<const Mesh> rMesh,
                                  BVH::View bvh,
                                  std::span<const CellData> contiguousCellData,
                                  Ref<const utils::ArgParse::Results> rArguments) {
    BoundaryMesh boundary;

    // Parse user input.
    const Scalar radius                     = rArguments.get<double>("boundary-radius");
    const std::size_t resolution            = rArguments.get<std::size_t>("boundary-resolution");
    const std::size_t minBoundaryTreeDepth  = rArguments.get<std::size_t>("min-boundary-tree-depth");
    const std::size_t maxBoundaryTreeDepth  = rArguments.get<std::size_t>("max-boundary-tree-depth");

    // Sanity checks.
    if (resolution < 3) CIE_THROW(Exception, "boundary resolution must be 3 or greater")

    // Generate vertices for the non-conforming boundary mesh.
    using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;
    DynamicArray<Point> corners;
    for (unsigned iSegment=0u; iSegment<resolution + 1; ++iSegment) {
        const Scalar arcParameter = iSegment * 2 * std::numbers::pi / resolution;
        corners.push_back(Point {
            radius * std::cos(arcParameter) + 0.5,
            radius * std::sin(arcParameter) + 0.5
        });
    } // for iSegment in range(resolution)

    // Generate edges for the non-conforming boundary mesh,
    // and cut them by the cells they're located in.
    for (unsigned iCorner=0u; iCorner<corners.size(); ++iCorner) {
        const std::array<Point,2> transformed {
            corners[iCorner],
            corners[(iCorner + 1) % corners.size()]};
        maths::AffineEmbedding<Scalar,1u,Dimension> transform(transformed);

        DynamicArray<ParametricBoundarySegment> boundarySegments;
        partitionBoundaryCell(
            transform,
            bvh,
            contiguousCellData,
            minBoundaryTreeDepth,
            maxBoundaryTreeDepth,
            boundarySegments);

        for (Ref<const ParametricBoundarySegment> rSegment : boundarySegments) {
            std::array<Point,2> segmentEndPoints;
            transform.evaluate(
                {&rSegment.segmentBegin, 1},
                segmentEndPoints.front());
            transform.evaluate(
                {&rSegment.segmentEnd, 1},
                segmentEndPoints.back());
            boundary.insert(BoundaryMesh::Vertex(
                boundary.vertices().size(),
                {},
                BoundaryCellData(
                    maths::AffineEmbedding<Scalar,1u,Dimension>(segmentEndPoints),
                    contiguousCellData[rSegment.iCell]
                )));
        } // for rSegment in boundarySegments
    } // for iCorner in range(1, corners.size())

    // Generate quadrature points.
    CIE_BEGIN_EXCEPTION_TRACING
        GaussLegendreQuadrature<Scalar> quadrature(boundaryIntegrationOrder);
        DynamicArray<QuadraturePoint<1,Scalar>> points;
        points.reserve(quadrature.numberOfNodes());
        for (std::size_t iPoint=0ul; iPoint<quadrature.numberOfNodes(); ++iPoint) {
            points.emplace_back(
                quadrature.nodes()[iPoint],
                quadrature.weights()[iPoint]);
        }
        boundary.data() = BoundaryMeshData(points, rMesh.data());
    CIE_END_EXCEPTION_TRACING

    return boundary;
}


using BoundarySegment = StaticArray<Scalar,2*Dimension>;


[[nodiscard]] DynamicArray<BoundarySegment>
imposeBoundaryConditions(Ref<Mesh> rMesh,
                         Ref<const Assembler> rAssembler,
                         BVH::View bvh,
                         std::span<const CellData> contiguousCellData,
                         CSRWrapper lhs,
                         std::span<Scalar> rhs,
                         Ref<const utils::ArgParse::Results> rArguments) {
    DynamicArray<BoundarySegment> boundarySegments;
    auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

    // Load the boundary mesh.
    const auto boundary = generateBoundaryMesh(
        rMesh,
        bvh,
        contiguousCellData,
        rArguments);

    const Quadrature<Scalar,1> lineQuadrature((GaussLegendreQuadrature<Scalar>(
        boundaryIntegrationOrder,
        utils::Comparison<Scalar>(1e-14, 1e-14),
        /*maxNewtonIterations=*/200ul)));
    DynamicArray<Scalar> quadratureBuffer;
    DynamicArray<Scalar> integrandBuffer;
    DirichletBoundary dirichletBoundary;
    const Scalar penaltyFactor = rArguments.get<double>("penalty-factor");
    const Scalar minBoundarySegmentNorm = rArguments.get<double>("min-boundary-segment-norm");

    for (const auto& rBoundaryCell : boundary.vertices()) {
        // Transform sample points to the global coordinate space.
        const std::array<ParametricCoordinate<Scalar>,1> segmentBegin {-1.0}, segmentEnd{1.0};
        std::array<
            std::array<
                PhysicalCoordinate<Scalar>,
                Dimension>,
            2
        > physicalCorners;
        rBoundaryCell.data().transform(
            segmentBegin,
            physicalCorners.front());
        rBoundaryCell.data().transform(
            segmentEnd,
            physicalCorners.back());

        // Discard the segment if it's too short.
        const Scalar segmentNorm = std::pow(physicalCorners.back()[0] - physicalCorners.front()[0], static_cast<Scalar>(2))
                                 + std::pow(physicalCorners.back()[1] - physicalCorners.front()[1], static_cast<Scalar>(2));

        if (minBoundarySegmentNorm < segmentNorm) {
            const auto ansatzSpace = rMesh.data().ansatzSpace();
            const maths::AffineEmbedding<Scalar,1,Dimension> segmentTransform(physicalCorners);

            auto integrand = makeTransformedIntegrand(
                makeDirichletPenaltyIntegrand(
                    dirichletBoundary,
                    penaltyFactor,
                    ansatzSpace,
                    segmentTransform,
                    rBoundaryCell.data().getEmbeddingCell()),
                segmentTransform.makeDerivative());
            integrandBuffer.resize(integrand.getMinBufferSize());
            integrand.setBuffer(integrandBuffer);
            quadratureBuffer.resize(integrand.size());
            lineQuadrature.evaluate(integrand, quadratureBuffer);

            const std::size_t lhsEntryCount = std::pow(ansatzSpace.size(), Dimension);
            rAssembler.addContribution<tags::Serial>(
                std::span<const Scalar>(quadratureBuffer.data(), lhsEntryCount),
                rBoundaryCell.data().getEmbeddingCell().id(),
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries);
            rAssembler.addContribution<tags::Serial>(
                std::span<const Scalar>(
                    quadratureBuffer.data() + lhsEntryCount,
                    quadratureBuffer.data() + quadratureBuffer.size()),
                rBoundaryCell.data().getEmbeddingCell().id(),
                std::span<Scalar>(rhs));

            // Log debug and output info.
            BoundarySegment segment;
            std::copy_n(physicalCorners[0].data(), Dimension, segment.data());
            std::copy_n(physicalCorners[1].data(), Dimension, segment.data() + Dimension);
            boundarySegments.push_back(segment);
        } // if minBoundarySegmentNorm < segmentNorm
    } // for rBoundaryCell in boundary.vertices()

    return boundarySegments;
}


} // namespace cie::fem
