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
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- STL Includes ---
#include <numbers> // std::numbers::pi
#include <regex>


namespace cie::fem {


/// @brief Data structure unique to the triangulated, immersed boundary cells.
class BoundaryCellData : public CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension> {
public:
    using Base = CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension>;

    using Base::Base;

    BoundaryCellData(unsigned id,
                     RightRef<typename Base::SpatialTransform> rEmbedding,
                     Ref<const std::array<Scalar,2>> state,
                     Ref<const CellData> rCell) noexcept
        : Base(
            VertexID(id),
            OrientedAxes<1>(),
            std::move(rEmbedding)),
         _pCell(&rCell),
         _state(state)
    {}

    Ref<const CellData> getEmbeddingCell() const noexcept {
        return *_pCell;
    }

    constexpr std::span<const Scalar,2> state() const noexcept {
        return std::span<const Scalar,2>(_state.data(), 2);
    }

private:
    Ptr<const CellData> _pCell;

    std::array<Scalar,2> _state;
}; // class BoundaryCellData


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
    void,
    BoundaryMeshData
>;


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


struct DirichletBoundary : public maths::ExpressionTraits<Scalar> {
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;

    constexpr DirichletBoundary() noexcept = default;

    constexpr DirichletBoundary(Ref<const std::span<const Scalar,2>> rState) noexcept
        : _state()
    {
        std::copy_n(
            rState.data(),
            rState.size(),
            _state.data());
    }

    void evaluate([[maybe_unused]] ConstSpan position, Span state) const noexcept {
        state.front() = _state.front() * (1.0 - 0.5 * (position.front() + 1.0)) + _state.back() * (0.5 * (position.front() + 1.0));
    }

    unsigned size() const noexcept {
        return 1;
    }

private:
    std::array<Scalar,2> _state;
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


using BoundarySegment = StaticArray<Scalar,2*Dimension+2>;


BoundaryMesh generateBoundaryMesh(BVH::View bvh,
                                  std::span<const CellData> contiguousCellData,
                                  Ref<const utils::ArgParse::Results> rArguments,
                                  Ref<DynamicArray<BoundarySegment>> rBoundarySegments) {
    BoundaryMesh boundary;

    // Parse user input.
    const std::size_t minBoundaryTreeDepth      = rArguments.get<std::size_t>("min-boundary-tree-depth");
    const std::size_t maxBoundaryTreeDepth      = rArguments.get<std::size_t>("max-boundary-tree-depth");
    const Scalar minBoundarySegmentNorm         = rArguments.get<double>("min-boundary-segment-norm");
    const std::filesystem::path boundaryFile    = rArguments.get<std::filesystem::path>("boundary-file-path");

    std::vector<BoundarySegment> segments;
    {
        std::ifstream file(boundaryFile);
        const std::string floatingPointRegex(R"(-?(?:(?:(?:[1-9][0-9]*)(?:\.[0-9]*)?)|(?:0(?:\.[0-9]*)?))(?:[eE][\+-]?[0-9]+)?)");
        const std::regex pattern(std::format(
            "^({}),({}),({}),({}),({}),({}).*",
            floatingPointRegex, floatingPointRegex, floatingPointRegex,
            floatingPointRegex, floatingPointRegex, floatingPointRegex));
        std::string line, component;
        while (std::getline(file, line)) {
            std::match_results<std::string::iterator> match;
            if (std::regex_match(line.begin(), line.end(), match, pattern)) {
                CIE_CHECK(match.size() == 6 + 1, "invalid line in boundary file: '" << line << "'" << "(" << match.size() << " matches)")
                BoundarySegment segment;
                std::transform(
                    match.begin() + 1,
                    match.end(),
                    segment.begin(),
                    [] (const auto& rSubMatch) -> Scalar {
                        const std::string& rString = rSubMatch.str();
                        return static_cast<Scalar>(std::stold(rString));
                    });
                segments.push_back(segment);
            } // if regex_match
        } // while getline
    }

    //// Generate vertices for the non-conforming boundary mesh.
    //using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;
    //DynamicArray<Point> corners;
    //for (unsigned iSegment=0u; iSegment<resolution + 1; ++iSegment) {
    //    const Scalar arcParameter = iSegment * 2 * std::numbers::pi / resolution;
    //    corners.push_back(Point {
    //        radius * std::cos(arcParameter) + Scalar(0.5),
    //        radius * std::sin(arcParameter) + Scalar(0.5)
    //    });
    //} // for iSegment in range(resolution)

    // Generate edges for the non-conforming boundary mesh,
    // and cut them by the cells they're located in.
    using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;
    for (Ref<const BoundarySegment> rPhysicalSegment : segments) {
        const std::array<Point,2> transformed {
            Point {rPhysicalSegment[0], rPhysicalSegment[1]},
            Point {rPhysicalSegment[2], rPhysicalSegment[3]}};
        maths::AffineEmbedding<Scalar,1u,Dimension> transform(transformed);

        DynamicArray<ParametricBoundarySegment> boundarySegments;
        partitionBoundaryCell(
            transform,
            bvh,
            contiguousCellData,
            minBoundaryTreeDepth,
            maxBoundaryTreeDepth,
            boundarySegments);

        for (Ref<const ParametricBoundarySegment> rParametricSegment : boundarySegments) {
            std::array<Point,2> segmentEndPoints;
            transform.evaluate(
                {&rParametricSegment.segmentBegin, 1},
                segmentEndPoints.front());
            transform.evaluate(
                {&rParametricSegment.segmentEnd, 1},
                segmentEndPoints.back());

            const Scalar segmentNorm = std::pow(segmentEndPoints.back()[0] - segmentEndPoints.front()[0], static_cast<Scalar>(2))
                                     + std::pow(segmentEndPoints.back()[1] - segmentEndPoints.front()[1], static_cast<Scalar>(2));

            if (minBoundarySegmentNorm < segmentNorm) {
                const auto id = boundary.vertices().size();
                const std::array<Scalar,2> state {
                    rPhysicalSegment[4] * (1.0 - 0.5 * (rParametricSegment.segmentBegin + 1.0)) + rPhysicalSegment[5] * (0.5 * (rParametricSegment.segmentBegin + 1.0)),
                    rPhysicalSegment[4] * (1.0 - 0.5 * (rParametricSegment.segmentEnd + 1.0)) + rPhysicalSegment[5] * (0.5 * (rParametricSegment.segmentEnd + 1.0))};
                boundary.insert(BoundaryMesh::Vertex(
                    id,
                    {},
                    BoundaryCellData(
                        id,
                        maths::AffineEmbedding<Scalar,1u,Dimension>(segmentEndPoints),
                        state,
                        contiguousCellData[rParametricSegment.iCell]
                    )));
                rBoundarySegments.push_back({
                    segmentEndPoints.front().front(),
                    segmentEndPoints.front().back(),
                    segmentEndPoints.back().front(),
                    segmentEndPoints.back().back(),
                    state.front(),
                    state.back()});
            }
        } // for rParametricSegment in boundarySegments
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
        boundary.data() = BoundaryMeshData(points);
    CIE_END_EXCEPTION_TRACING

    return boundary;
}


[[nodiscard]] DynamicArray<BoundarySegment>
imposeBoundaryConditions(Ref<Mesh> rMesh,
                         Ref<const Assembler> rAssembler,
                         BVH::View bvh,
                         std::span<const CellData> contiguousCellData,
                         CSRWrapper lhs,
                         std::span<Scalar> rhs,
                         Ref<const utils::ArgParse::Results> rArguments) {
    DynamicArray<BoundarySegment> boundarySegments;

    CIE_BEGIN_EXCEPTION_TRACING
    auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

    // Load the boundary mesh.
    const auto boundary = generateBoundaryMesh(
        bvh,
        contiguousCellData,
        rArguments,
        boundarySegments);

    const auto quadratureRuleFactory = [&boundary] (Ref<const BoundaryMesh::Vertex::Data>) {
        return boundary.data().makeQuadratureRule();};

    const Scalar penaltyFactor = rArguments.get<double>("penalty-factor");
    DynamicArray<Scalar> integrandBuffer;
    const auto integrandFactory = [&rMesh, penaltyFactor, &integrandBuffer] (Ref<const BoundaryMesh::Vertex::Data> rBoundaryCell) {
        auto integrand = makeTransformedIntegrand(
            makeDirichletPenaltyIntegrand(
                DirichletBoundary(rBoundaryCell.state()),
                penaltyFactor,
                rMesh.data().ansatzSpace(),
                rBoundaryCell.makeSpatialTransform(),
                rBoundaryCell.getEmbeddingCell()),
            rBoundaryCell.makeJacobian());
        integrandBuffer.resize(integrand.getMinBufferSize());
        integrand.setBuffer(integrandBuffer);
        return integrand;
    }; // integrandFactory

    const auto integralSink = [&lhs, &rhs, &rAssembler, &boundary] (std::span<const VertexID> cellIDs,
                                                                    std::span<const Scalar> results) {
        constexpr std::size_t lhsEntryCount = intPow(Ansatz::size(), Dimension);
        constexpr std::size_t rhsEntryCount = Ansatz::size();
        for (std::size_t iCell=0ul; iCell<cellIDs.size(); ++iCell) {
            Ptr<const Scalar> pResultsBegin = results.data() + iCell * (lhsEntryCount +rhsEntryCount);
            const VertexID cellID = boundary.find(cellIDs[iCell]).value().data().getEmbeddingCell().id();
            rAssembler.addContribution<tags::Serial>(
                std::span<const Scalar>(pResultsBegin, lhsEntryCount),
                cellID,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries);
            rAssembler.addContribution<tags::Serial>(
                std::span<const Scalar>(pResultsBegin + lhsEntryCount, rhsEntryCount),
                cellID,
                std::span<Scalar>(rhs));
        } // for iCell in range(cellIDs.size())
    };

    using Integrand = TransformedIntegrand<
        DirichletPenaltyIntegrand<
            DirichletBoundary,
            Ansatz,
            maths::AffineEmbedding<Scalar,1u,Dimension>,
            CellData>,
        maths::AffineEmbedding<Scalar,1u,Dimension>::Derivative>;
    const IntegrandProcessor<1,Integrand>::Properties executionProperties{
        .integrandBatchSize = rArguments.get<std::size_t>("integrand-batch-size"),
        .integrandsPerItem = {},
        .verbosity = 3};
    auto pProcessor = std::make_unique<IntegrandProcessor<1,Integrand>>();
    const auto& rBoundaryCells = boundary.vertices();
    pProcessor->process(
        rBoundaryCells.begin(),
        rBoundaryCells.end(),
        quadratureRuleFactory,
        integrandFactory,
        integralSink,
        executionProperties);

    return boundarySegments;
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
