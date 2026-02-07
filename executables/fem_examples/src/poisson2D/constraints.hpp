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

    BoundaryCellData(RightRef<typename Base::SpatialTransform> rTransform) noexcept
        : Base(
            VertexID(),
            Base::AnsatzSpaceID(),
            OrientedAxes<1>(),
            std::move(rTransform))
    {}
}; // class BoundaryCellData
//using BoundaryCellData = maths::AffineEmbedding<Scalar,1u,Dimension>;


/// @brief Data structure unqiue to the triangulated, immersed boundary cell corners.
using BoundaryCornerData = Scalar;


class BoundaryMeshData {
public:
    BoundaryMeshData() noexcept = default;

    BoundaryMeshData(std::span<const QuadraturePoint<1,Scalar>> quadraturePointSet)
        : _quadraturePointSet(quadraturePointSet.begin(), quadraturePointSet.end())
    {}

    CachedQuadraturePointFactory<1,Scalar> makeQuadratureRule() const {
        return CachedQuadraturePointFactory<1,Scalar>(
            std::span<const QuadraturePoint<1,Scalar>>(
                _quadraturePointSet.data(),
                _quadraturePointSet.size()));
    }

private:
    std::vector<QuadraturePoint<1,Scalar>> _quadraturePointSet;
}; // class BoundaryMeshData


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    BoundaryCornerData,
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


BoundaryMesh generateBoundaryMesh(const Scalar radius, const unsigned resolution) {
    using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;

    if (resolution < 3) CIE_THROW(Exception, "boundary resolution must be 3 or greater")

    BoundaryMesh boundary;

    DynamicArray<Point> corners;
    for (unsigned iSegment=0u; iSegment<resolution + 1; ++iSegment) {
        const Scalar arcParameter = iSegment * 2 * std::numbers::pi / resolution;
        corners.push_back(Point {
            radius * std::cos(arcParameter) + 0.5,
            radius * std::sin(arcParameter) + 0.5
        });
    } // for iSegment in range(resolution)

    for (unsigned iCorner=0u; iCorner<corners.size(); ++iCorner) {
        const std::array<Point,2> transformed {
            corners[iCorner],
            corners[(iCorner + 1) % corners.size()]};

        boundary.insert(BoundaryMesh::Vertex(
            boundary.vertices().size(),
            {},
            {maths::AffineEmbedding<Scalar,1u,Dimension>(transformed)}
        ));
    } // for iCorner in range(1, corners.size())

    for (unsigned iCorner=0; iCorner<corners.size(); ++iCorner) {
        boundary.insert(BoundaryMesh::Edge(
            iCorner,
            {iCorner, (iCorner + 1) % boundary.vertices().size()},
            static_cast<Scalar>(iCorner)
        ));
    } // for iCorner in range(corners.size())

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
        boundary.data() = BoundaryMeshData(points);
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
    auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

    DynamicArray<BoundarySegment> boundarySegments; // {p0x, p0y, p1x, p1y}
    const unsigned boundaryResolution = rArguments.get<std::size_t>("boundary-resolution");

    // Load the boundary mesh.
    const Scalar boundaryRadius = rArguments.get<double>("boundary-radius");
    const auto boundary = generateBoundaryMesh(boundaryRadius, boundaryResolution);

    using TreePrimitive = geo::Cube<1,Scalar>;
    using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

    const Quadrature<Scalar,1> lineQuadrature((GaussLegendreQuadrature<Scalar>(
        boundaryIntegrationOrder,
        utils::Comparison<Scalar>(1e-14, 1e-14),
        /*maxNewtonIterations=*/200ul)));
    DynamicArray<Scalar> quadratureBuffer;
    DynamicArray<Scalar> integrandBuffer;
    DirichletBoundary dirichletBoundary;
    const Scalar penaltyFactor = rArguments.get<double>("penalty-factor");
    const unsigned minBoundaryTreeDepth = rArguments.get<std::size_t>("min-boundary-tree-depth");
    const unsigned maxBoundaryTreeDepth = rArguments.get<std::size_t>("max-boundary-tree-depth");
    const Scalar minBoundarySegmentNorm = rArguments.get<double>("min-boundary-segment-norm");

    for (const auto& rBoundaryCell : boundary.vertices()) {
        Tree tree(Tree::Point {-1.0}, 2.0);
        std::vector<std::tuple<
            std::size_t,    // cell data index
            Scalar,         // segment begin
            Scalar          // segment end
        >> homogeneousSegments;

        // Define a functor detecting cell boundaries.
        const auto boundaryVisitor = [bvh, &tree, &rBoundaryCell,
                                      &contiguousCellData,
                                      &homogeneousSegments,
                                      minBoundaryTreeDepth, maxBoundaryTreeDepth] (
            Ref<const Tree::Node> rNode,
            unsigned level) -> bool {

            if (level < minBoundaryTreeDepth) return true;
            if (maxBoundaryTreeDepth < level) return false;

            Scalar base, edgeLength;
            tree.getNodeGeometry(rNode, &base, &edgeLength);

            // Define sample points in the boundary cell's parametric space.
            const std::array<ParametricCoordinate<Scalar>,1>
                parametricBase {base},
                parametricOpposite {base + edgeLength};

            // Transform sample points to the global coordinate space.
            std::array<PhysicalCoordinate<Scalar>,Dimension> physicalBase, physicalOpposite;
            rBoundaryCell.data().transform(parametricBase, physicalBase);
            rBoundaryCell.data().transform(parametricOpposite, physicalOpposite);

            // Check whether the two endpoints are in different cells.
            const auto iBaseCell = bvh.find(
                std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                    physicalBase.data()),
                    Dimension),
                contiguousCellData);
            const auto iOppositeCell = bvh.find(
                std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                    physicalOpposite.data()),
                    Dimension),
                contiguousCellData);

            // Integrate if both endpoints lie in the same cell.
            if (iBaseCell != contiguousCellData.size() && iOppositeCell != contiguousCellData.size() && iBaseCell == iOppositeCell) {
                homogeneousSegments.emplace_back(
                    iBaseCell,
                    base,
                    base + edgeLength);
            } // if both endpoints lie in the same cell

            return iBaseCell != iOppositeCell &&
                   (iBaseCell != contiguousCellData.size() || iOppositeCell != contiguousCellData.size());
        }; // boundaryVisitor

        // Construct a binary tree that detects intersections between
        // Cell boundaries and the current boundary cell.
        tree.scan(boundaryVisitor);

        std::sort(
            homogeneousSegments.begin(),
            homogeneousSegments.end(),
            [](const auto& rLeft, const auto& rRight) {
                if (std::get<0>(rLeft) < std::get<0>(rRight)) return true;
                else if (std::get<0>(rRight) < std::get<0>(rLeft)) return false;
                else return std::get<1>(rLeft) < std::get<1>(rRight);
            });

        auto itBegin = homogeneousSegments.begin();
        while (itBegin != homogeneousSegments.end()) {
            const std::size_t iCellData = std::get<0>(*itBegin);
            const Scalar segmentBegin = std::get<1>(*itBegin);

            auto itEnd = std::upper_bound(
                itBegin,
                homogeneousSegments.end(),
                iCellData,
                [](std::size_t iCellData, const auto& rItem) {
                    return iCellData < std::get<0>(rItem);
                });
            const Scalar segmentEnd = itBegin == itEnd ? std::get<2>(*itEnd) : std::get<2>(*(itEnd - 1));

            // Transform sample points to the global coordinate space.
            std::array<
                std::array<
                    PhysicalCoordinate<Scalar>,
                    Dimension>,
                2
            > physicalCorners;
            rBoundaryCell.data().transform(
                Kernel<1,Scalar>::cast<ParametricCoordinate<Scalar>>(std::span<const Scalar,1>(&segmentBegin, 1)),
                physicalCorners.front());
            rBoundaryCell.data().transform(
                Kernel<1,Scalar>::cast<ParametricCoordinate<Scalar>>(std::span<const Scalar,1>(&segmentEnd, 1)),
                physicalCorners.back());

            // Discard the segment if it's too short.
            const Scalar segmentNorm =   std::pow(physicalCorners.back()[0] - physicalCorners.front()[0], static_cast<Scalar>(2))
                                       + std::pow(physicalCorners.back()[1] - physicalCorners.front()[1], static_cast<Scalar>(2));

            if (minBoundarySegmentNorm < segmentNorm) {
                Ref<const CellData> rCell = contiguousCellData[iCellData];
                const auto ansatzSpace = rMesh.data().ansatzSpace();
                const maths::AffineEmbedding<Scalar,1,Dimension> segmentTransform(physicalCorners);

                auto integrand = makeTransformedIntegrand(
                    makeDirichletPenaltyIntegrand(
                        dirichletBoundary,
                        penaltyFactor,
                        ansatzSpace,
                        segmentTransform,
                        rCell),
                    segmentTransform.makeDerivative());
                integrandBuffer.resize(integrand.getMinBufferSize());
                integrand.setBuffer(integrandBuffer);
                quadratureBuffer.resize(integrand.size());
                lineQuadrature.evaluate(integrand, quadratureBuffer);

                const std::size_t lhsEntryCount = std::pow(ansatzSpace.size(), Dimension);
                rAssembler.addContribution<tags::Serial>(
                    std::span<const Scalar>(quadratureBuffer.data(), lhsEntryCount),
                    rCell.id(),
                    lhs.rowExtents,
                    lhs.columnIndices,
                    lhs.entries);
                rAssembler.addContribution<tags::Serial>(
                    std::span<const Scalar>(
                        quadratureBuffer.data() + lhsEntryCount,
                        quadratureBuffer.data() + quadratureBuffer.size()),
                    rCell.id(),
                    std::span<Scalar>(rhs));

                // Log debug and output info.
                BoundarySegment segment;
                std::copy_n(physicalCorners[0].data(), Dimension, segment.data());
                std::copy_n(physicalCorners[1].data(), Dimension, segment.data() + Dimension);
                boundarySegments.push_back(segment);
            }

            itBegin = itEnd;
        } // while itBegin != homogeneousSegments.end()

    } // for rBoundaryCell in boundary.vertices()

    return boundarySegments;
}


} // namespace cie::fem
