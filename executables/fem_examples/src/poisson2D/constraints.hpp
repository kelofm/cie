#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/numeric/inc/Quadrature.hpp"

// --- STL Includes ---
#include <numbers> // std::numbers::pi


namespace cie::fem {


/// @brief Data structure unique to the triangulated, immersed boundary cells.
class BoundaryCellData : public CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,2u>,void,2u> {
public:
    using Base = CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,2u>,void,2u>;

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
    using DomainData = std::uint8_t;

private:
}; // class BoundaryMeshData


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    BoundaryCornerData
>;


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


struct DirichletBoundary : public maths::ExpressionTraits<Scalar> {
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;

    void evaluate(ConstSpan position, Span state) const noexcept {
        state[0] = position[0] + position[1];
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


BoundaryMesh generateBoundaryMesh(const unsigned resolution) {
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


[[nodiscard]] DynamicArray<StaticArray<Scalar,2*Dimension+1>>
imposeBoundaryConditions(Ref<Mesh> rMesh,
                         Ref<const Assembler> rAssembler,
                         BVH::View bvh,
                         std::span<const CellData> contiguousCellData,
                         CSRWrapper lhs,
                         std::span<Scalar> rhs,
                         Ref<const utils::ArgParse::Results> rArguments) {
    auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

    DynamicArray<StaticArray<Scalar,2*Dimension+1>> boundarySegments; // {p0x, p0y, p1x, p1y, level}
    const unsigned boundaryResolution = rArguments.get<std::size_t>("boundary-resolution");

    // Load the boundary mesh.
    const auto boundary = generateBoundaryMesh(boundaryResolution);

    using TreePrimitive = geo::Cube<1,Scalar>;
    using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

    const Quadrature<Scalar,1> lineQuadrature((GaussLegendreQuadrature<Scalar>(
        boundaryIntegrationOrder,
        utils::Comparison<Scalar>(1e-14, 1e-14),
        /*maxNewtonIterations=*/200ul)));
    DynamicArray<Scalar> quadratureBuffer(  std::pow(rMesh.data().ansatzSpace().size(), Dimension)
                                          + rMesh.data().ansatzSpace().size());
    DynamicArray<Scalar> integrandBuffer;

    DirichletBoundary dirichletBoundary;
    const Scalar penaltyFactor = rArguments.get<double>("penalty-factor");

    for (const auto& rBoundaryCell : boundary.vertices()) {
        Tree tree(Tree::Point {-1.0}, 2.0);

        // Define a functor detecting cell boundaries.
        const auto boundaryVisitor = [bvh, &tree, &rBoundaryCell,
                                      &lineQuadrature, &integrandBuffer, &quadratureBuffer, &rMesh,
                                      lhs, &rAssembler,
                                      &rhs, &dirichletBoundary,
                                      &boundarySegments, &contiguousCellData, penaltyFactor] (
            Ref<const Tree::Node> rNode,
            unsigned level) -> bool {

            if (level < minBoundaryTreeDepth) return true;
            if (maxBoundaryTreeDepth < level) return false;

            Scalar base, edgeLength;
            tree.getNodeGeometry(rNode, &base, &edgeLength);

            // Define sample points in the boundary cell's local space.
            const std::array<Kernel<1, Scalar>::LocalCoordinate,1>
                localBase {base},
                localOpposite {base + edgeLength};

            // Transform sample points to the global coordinate space.
            std::array<Kernel<Dimension, Scalar>::GlobalCoordinate,Dimension> globalBase, globalOpposite;
            rBoundaryCell.data().transform(localBase, globalBase);
            rBoundaryCell.data().transform(localOpposite, globalOpposite);

            // Check whether the two endpoints are in different cells.
            const auto iBaseCell = bvh.find(
                std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(globalBase.data()), Dimension),
                contiguousCellData);
            const auto iOppositeCell = bvh.find(
                std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(globalOpposite.data()), Dimension),
                contiguousCellData);

            // Integrate if both endpoints lie in the same cell.
            if (iBaseCell != contiguousCellData.size() && iOppositeCell != contiguousCellData.size() && iBaseCell == iOppositeCell) {
                const Scalar segmentNorm =   std::pow(globalOpposite[0] - globalBase[0], static_cast<Scalar>(2))
                                           + std::pow(globalOpposite[1] - globalBase[1], static_cast<Scalar>(2));

                if (minBoundarySegmentNorm < segmentNorm) {
                    Ref<const CellData> rCell = contiguousCellData[iBaseCell];
                    const auto ansatzSpace = rMesh.data().ansatzSpace();

                    StaticArray<maths::AffineEmbedding<Scalar,1,Dimension>::OutPoint,2> globalCorners;
                    globalCorners[0][0] = globalBase[0];
                    globalCorners[0][1] = globalBase[1];
                    globalCorners[1][0] = globalOpposite[0];
                    globalCorners[1][1] = globalOpposite[1];
                    const maths::AffineEmbedding<Scalar,1,Dimension> segmentTransform(globalCorners);

                    auto integrand = makeTransformedIntegrand(
                        makeDirichletPenaltyIntegrand(
                            dirichletBoundary,
                            penaltyFactor,
                            ansatzSpace,
                            segmentTransform,
                            std::span<Scalar>(integrandBuffer)),
                        segmentTransform.makeDerivative());
                    integrandBuffer.resize(integrand.getMinBufferSize());
                    integrand.setBuffer(integrandBuffer);
                    lineQuadrature.evaluate(integrand, quadratureBuffer);

                    const std::size_t lhsEntryCount = std::pow(ansatzSpace.size(), Dimension);
                    rAssembler.addContribution(
                        std::span<const Scalar>(quadratureBuffer.data(), lhsEntryCount),
                        rCell.id(),
                        lhs.rowExtents,
                        lhs.columnIndices,
                        lhs.entries);
                    rAssembler.addContribution(
                        std::span<const Scalar>(quadratureBuffer.data() + lhsEntryCount, quadratureBuffer.size() - lhsEntryCount),
                        rCell.id(),
                        std::span<Scalar>(rhs));

                    // Log debug and output info.
                    decltype(boundarySegments)::value_type segment;
                    std::copy_n(globalCorners[0].data(), Dimension, segment.data());
                    std::copy_n(globalCorners[1].data(), Dimension, segment.data() + Dimension);
                    segment.back() = level;
                    boundarySegments.push_back(segment);
                }
            } // if both endpoints lie in the same cell

            return iBaseCell != iOppositeCell &&
                   (iBaseCell != contiguousCellData.size() || iOppositeCell != contiguousCellData.size());
        }; // boundaryVisitor

        // Construct a binary tree that detects intersections between
        // Cell boundaries and the current boundary cell.
        tree.scan(boundaryVisitor);
    } // for rBoundaryCell in boundary.vertices()

    return boundarySegments;
}


} // namespace cie::fem
