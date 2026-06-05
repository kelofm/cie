// --- Internal Includes ---
#include "poisson2D/constraints.hpp"

// --- FEM Includes ---
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- GEO Includes ---
#include "packages/trees/inc/ContiguousSpaceTree.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/io/inc/json.hpp"

// --- STL Includes ---
#include <regex>


namespace cie::fem {


BoundaryCell::BoundaryCell(
    unsigned id,
    std::size_t ansatzID,
    RightRef<typename Base::SpatialTransform> rEmbedding,
    Ref<const std::array<Scalar,2>> state,
    Ref<const Cell> rCell) noexcept
        : Base(
            id,
            ansatzID,
            OrientedAxes<1>(),
            std::move(rEmbedding)),
          _pCell(&rCell),
          _state(state)
{}


BoundaryMeshData::BoundaryMeshData(
    std::span<const QuadraturePoint<1,Scalar>> quadraturePointSet,
    RightRef<std::vector<BoundaryCell>> rBoundaryCells)
        :   _quadraturePointSet(quadraturePointSet.begin(), quadraturePointSet.end()),
            _cells(std::move(rBoundaryCells))
{}


CachedQuadraturePointFactory<1,Scalar> BoundaryMeshData::makeQuadratureRule() const {
    return CachedQuadraturePointFactory<1,Scalar>(
        std::span<const QuadraturePoint<1,Scalar>>(
            _quadraturePointSet.data(),
            _quadraturePointSet.size()));
}


struct DirichletBoundary : public maths::ExpressionTraits<Scalar> {
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;
    using maths::ExpressionTraits<Scalar>::BufferSpan;

    constexpr DirichletBoundary() noexcept = default;

    constexpr DirichletBoundary(Ref<const std::span<const Scalar,2>> rState) noexcept
        : _state() {
            std::copy_n(
                rState.data(),
                rState.size(),
                _state.data());
    }

    void evaluate(
        ConstSpan position,
        Span state,
        BufferSpan) const noexcept {
            state.front() = _state.front() * (1.0 - 0.5 * (position.front() + 1.0)) + _state.back() * (0.5 * (position.front() + 1.0));
    }

    unsigned size() const noexcept {
        return 1;
    }

    constexpr unsigned bufferSize() const noexcept {
        return 0u;
    }

private:
    std::array<Scalar,2> _state;
}; // class DirichletBoundary


/// @brief Definition of a boundary segment in parametric space.
struct ParametricBoundarySegment {
    std::size_t iCell;
    Scalar segmentBegin;
    Scalar segmentEnd;
}; // struct ParametricBoundarySegment


void partitionBoundaryCell(
    Ref<maths::AffineEmbedding<Scalar,1u,Dimension>> rTransform,
    BVH::View bvh,
    std::span<const Cell> cells,
    unsigned minBoundaryTreeDepth,
    unsigned maxBoundaryTreeDepth,
    Ref<DynamicArray<ParametricBoundarySegment>> rOutput) {
        rOutput.clear();
        CIE_BEGIN_EXCEPTION_TRACING
            using TreePrimitive = geo::Cube<1,Scalar>;
            using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;
            Tree tree(Tree::Point {-1.0}, 2.0);
            std::vector<Scalar> buffer;

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
                buffer.resize(rTransform.bufferSize());
                rTransform.evaluate(parametricBase, physicalBase, buffer);
                rTransform.evaluate(parametricOpposite, physicalOpposite, buffer);

                // Check whether the two endpoints are in different cells.
                const std::size_t iBaseCell = bvh.find(
                    std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                        physicalBase.data()),
                        Dimension),
                    cells);
                const std::size_t iOppositeCell = bvh.find(
                    std::span<const Scalar,Dimension>(reinterpret_cast<const Scalar*>(
                        physicalOpposite.data()),
                        Dimension),
                    cells);

                // Extend the output if both endpoints lie in the same cell.
                if (iBaseCell == iOppositeCell && iBaseCell != cells.size() && iOppositeCell != cells.size()) {
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
                    && (iBaseCell != cells.size() || iOppositeCell != cells.size());
            }; // visitor

            tree.scan(visitor);
        CIE_END_EXCEPTION_TRACING
}


std::vector<BoundarySegment> makeBoundary(Ref<const cie::io::JSONObject> configuration) {
    std::vector<BoundarySegment> output;
    CIE_BEGIN_EXCEPTION_TRACING
        const std::filesystem::path boundaryFile = configuration["file-path"].as<std::string>();
        std::cout << "reading " << boundaryFile << std::endl;
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
                CIE_CHECK(
                    match.size() == 6 + 1,
                    std::format(
                        "invalid lin in boundary file: '{}' ({} matches)",
                        line, match.size()))
                BoundarySegment segment;
                std::transform(
                    match.begin() + 1,
                    match.end(),
                    segment.begin(),
                    [] (const auto& rSubMatch) -> Scalar {
                        const std::string& rString = rSubMatch.str();
                        return static_cast<Scalar>(std::stold(rString));
                    });
                output.push_back(segment);
            } // if regex_match
        } // while getline
    CIE_END_EXCEPTION_TRACING

    return output;
}


BoundaryMesh generateBoundaryMesh(
    std::span<const BoundarySegment> tesselatedBoundary,
    BVH::View bvh,
    std::span<const Cell> cells,
    Ref<const cie::io::JSONObject> rConfiguration,
    Ref<DynamicArray<BoundarySegment>> rBoundarySegments) {
        BoundaryMesh boundary;
        std::vector<BoundaryCell> boundaryCells;

        // Parse user input.
        const std::size_t minTreeDepth  = rConfiguration["integration"]["min-tree-depth"].as<std::size_t>();
        const std::size_t maxTreeDepth  = rConfiguration["integration"]["max-tree-depth"].as<std::size_t>();
        const Scalar minSegmentNorm     = rConfiguration["integration"]["min-domain-norm"].as<std::size_t>();

        // Generate edges for the non-conforming boundary mesh,
        // and cut them by the cells they're located in.
        using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;
        std::vector<Scalar> buffer;

        for (Ref<const BoundarySegment> rPhysicalSegment : tesselatedBoundary) {
            const std::array<Point,2> transformed {
                Point {rPhysicalSegment[0], rPhysicalSegment[1]},
                Point {rPhysicalSegment[2], rPhysicalSegment[3]}};
            maths::AffineEmbedding<Scalar,1u,Dimension> transform(transformed);
            buffer.resize(transform.bufferSize());

            DynamicArray<ParametricBoundarySegment> boundarySegments;
            partitionBoundaryCell(
                transform,
                bvh,
                cells,
                minTreeDepth,
                maxTreeDepth,
                boundarySegments);

            for (Ref<const ParametricBoundarySegment> rParametricSegment : boundarySegments) {
                std::array<Point,2> segmentEndPoints;
                transform.evaluate(
                    {&rParametricSegment.segmentBegin, 1},
                    segmentEndPoints.front(),
                    buffer);
                transform.evaluate(
                    {&rParametricSegment.segmentEnd, 1},
                    segmentEndPoints.back(),
                    buffer);

                const Scalar segmentNorm = std::pow(segmentEndPoints.back()[0] - segmentEndPoints.front()[0], static_cast<Scalar>(2))
                                        + std::pow(segmentEndPoints.back()[1] - segmentEndPoints.front()[1], static_cast<Scalar>(2));

                if (minSegmentNorm < segmentNorm) {
                    const auto id = boundary.vertices().size();
                    const std::array<Scalar,2> state {
                        Scalar(rPhysicalSegment[4] * (1.0 - 0.5 * (rParametricSegment.segmentBegin + 1.0)) + rPhysicalSegment[5] * (0.5 * (rParametricSegment.segmentBegin + 1.0))),
                        Scalar(rPhysicalSegment[4] * (1.0 - 0.5 * (rParametricSegment.segmentEnd + 1.0)) + rPhysicalSegment[5] * (0.5 * (rParametricSegment.segmentEnd + 1.0)))};
                    boundary.insert(BoundaryMesh::Vertex(
                        id,
                        {}));
                    boundaryCells.emplace_back(
                        id,
                        0ul,
                        maths::AffineEmbedding<Scalar,1u,Dimension>(segmentEndPoints),
                        state,
                        cells[rParametricSegment.iCell]);
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
            GaussLegendreQuadrature<Scalar> quadrature(rConfiguration["integration"]["template"]["order"].as<std::size_t>());
            DynamicArray<QuadraturePoint<1,Scalar>> points;
            points.reserve(quadrature.numberOfNodes());
            for (std::size_t iPoint=0ul; iPoint<quadrature.numberOfNodes(); ++iPoint) {
                points.emplace_back(
                    quadrature.nodes()[iPoint],
                    quadrature.weights()[iPoint]);
            }
            boundary.data() = BoundaryMeshData(
                points,
                std::move(boundaryCells));
        CIE_END_EXCEPTION_TRACING

        return boundary;
}


DynamicArray<BoundarySegment>
integrateBoundaryConstraints(
    Ref<Mesh> rMesh,
    std::span<const BoundarySegment> tesselatedBoundary,
    Ref<const Assembler> rAssembler,
    BVH::View bvh,
    std::span<const Cell> cells,
    linalg::CSRView<Scalar,const int> lhs,
    Ref<DynamicArray<int>> rConstraintRowExtents,
    Ref<DynamicArray<int>> rConstraintColumnIndices,
    Ref<DynamicArray<Scalar>> rConstraintEntries,
    Ref<DynamicArray<Scalar>> rConstraintGaps,
    [[maybe_unused]] Ref<mp::ThreadPoolBase> rThreads,
    Ref<const cie::io::JSONObject> rConfiguration) {
        DynamicArray<BoundarySegment> boundarySegments;

        using Integrand = DirichletPenaltyIntegrand<
            DirichletBoundary,
            Ansatz,
            maths::AffineEmbedding<Scalar,1u,Dimension>,
            Cell>;

        CIE_BEGIN_EXCEPTION_TRACING
            auto logBlock = utils::LoggerSingleton::get().newBlock("weak boundary condition imposition");

            // Load the boundary mesh.
            const auto boundaryMesh = generateBoundaryMesh(
                tesselatedBoundary,
                bvh,
                cells,
                rConfiguration,
                boundarySegments);

            // Construct a new linear system containing only
            // contributions from constraints.
            std::vector<VertexID> constrainedCellIDs;
            {
                tsl::robin_set<VertexID> constrainedCellIDSet;
                for (Ref<const BoundaryCell> rBoundaryCell : boundaryMesh.data().cells())
                    constrainedCellIDSet.insert(rBoundaryCell.getEmbeddingCell().id());
                constrainedCellIDs.insert(
                    constrainedCellIDs.end(),
                    constrainedCellIDSet.begin(),
                    constrainedCellIDSet.end());
            }

            int rowCount, columnCount;
            rAssembler.makeCSRMatrix(
                rowCount,
                columnCount,
                rConstraintRowExtents,
                rConstraintColumnIndices,
                rConstraintEntries,
                constrainedCellIDs);

            // Extend the constraint system to match the unconstrained system's size.
            rConstraintRowExtents.resize(lhs.rowCount() + 1);
            std::fill(
                rConstraintRowExtents.begin() + rowCount,
                rConstraintRowExtents.end(),
                rConstraintRowExtents[rowCount]);
            rowCount = lhs.rowCount();
            columnCount = lhs.columnCount();

            linalg::CSRView<Scalar,int> constraintGradient(
                columnCount,
                rConstraintRowExtents,
                rConstraintColumnIndices,
                rConstraintEntries);
            rConstraintGaps.resize(rowCount);
            std::fill(
                constraintGradient.entries().begin(),
                constraintGradient.entries().end(),
                Scalar(0));
            std::fill(
                rConstraintGaps.begin(),
                rConstraintGaps.end(),
                Scalar(0));

            // Define callbacks for the integrand processor.
            const auto quadratureRuleFactory = [&boundaryMesh] (Ref<const BoundaryCell>) {
                return boundaryMesh.data().makeQuadratureRule();};

            const auto integrandFactory = [&rMesh] (Ref<const BoundaryCell> rBoundaryCell) -> Integrand {
                auto integrand = makeDirichletPenaltyIntegrand(
                        DirichletBoundary(rBoundaryCell.state()),
                        /*penaltyFactor=*/1,
                        rMesh.data().ansatz(rBoundaryCell.ansatzID()),
                        rBoundaryCell.makeSpatialTransform(),
                        rBoundaryCell.getEmbeddingCell());
                return integrand;
            }; // integrandFactory

            const unsigned integrandSize = integrandFactory(boundaryMesh.data().cells().front()).size();
            const unsigned constraintGradientContributionSize = intPow(Ansatz::size(), 2);
            const unsigned constraintGapContributionSize = Ansatz::size();
            const auto integralSink = [&] (
                std::span<const VertexID> cellIDs,
                std::span<const Scalar> results) {
                    for (std::size_t iCell=0ul; iCell<cellIDs.size(); ++iCell) {
                        // Find the results' offset.
                        Ptr<const Scalar> pResultsBegin = results.data() + iCell * integrandSize;
                        const VertexID embeddingCellID = boundaryMesh.data().cells()[iCell].getEmbeddingCell().id();

                        // Separate different integral components.
                        const std::span<const Scalar> constraintGradientContribution(
                            pResultsBegin,
                            constraintGradientContributionSize);
                        const std::span<const Scalar> constraintGapContribution(
                            pResultsBegin + constraintGradientContributionSize,
                            constraintGapContributionSize);

                        // Reduce contributions to the constraint system.
                        rAssembler.addContribution<tags::SMP,int>(
                            constraintGradientContribution,
                            embeddingCellID,
                            constraintGradient.rowExtents(),
                            constraintGradient.columnIndices(),
                            constraintGradient.entries());
                        rAssembler.addContribution<tags::SMP>(
                            constraintGapContribution,
                            embeddingCellID,
                            std::span<Scalar>(rConstraintGaps));
                    } // for iCell in range(cellIDs.size())
            }; // integralSink

            // Perform the integration.
            const IntegrandProcessor<1,Integrand>::Properties executionProperties{
                .integrandBatchSize = rConfiguration["integration"]["batch-size"].as<std::size_t>(),
                .verbosity = rConfiguration["integration"]["verbosity"].as<int>()};
            auto pProcessor = std::make_unique<IntegrandProcessor<1,Integrand>>();
            const auto& rBoundaryCells = boundaryMesh.data().cells();
            pProcessor->process(
                rBoundaryCells.begin(),
                rBoundaryCells.end(),
                quadratureRuleFactory,
                integrandFactory,
                integralSink,
                executionProperties);

//            // Shrink the constraint system.
//            {
//                std::vector<int> swapConstraintRowExtents(1, 0);
//                swapConstraintRowExtents.reserve(rConstraintColumnIndices.size());
//                std::vector<Scalar> swapConstraintGaps;
//                swapConstraintGaps.reserve(rConstraintGaps.size());
//
//                for (int iRow=0ul; iRow<lhs.rowCount(); ++iRow) {
//                    const int rowSize = rConstraintRowExtents[iRow + 1] - rConstraintRowExtents[iRow];
//                    if (rowSize) [[unlikely]] {
//                        swapConstraintRowExtents.push_back(rConstraintRowExtents[iRow + 1]);
//                        swapConstraintGaps.push_back(rConstraintGaps[iRow]);
//                    }
//                } // for iRow
//
//                rConstraintRowExtents = std::move(swapConstraintRowExtents);
//                rConstraintGaps = std::move(swapConstraintGaps);
//            }

            return boundarySegments;
        CIE_END_EXCEPTION_TRACING
}


BVH makeBoundingVolumeHierarchy(
    std::span<Cell> cells,
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths) {
        auto logBlock = utils::LoggerSingleton::get().newBlock("make BVH");

        constexpr int targetLeafWidth = 5;
        constexpr int maxTreeDepth = 5;

        geo::AABBoxNode<Cell> root;
        geo::AABBoxNode<Cell>::Point rootBase, rootLengths;
        std::transform(
            meshBase.begin(),
            meshBase.end(),
            meshLengths.begin(),
            rootBase.begin(),
            [] (Scalar base, Scalar length) -> Scalar {
                return base - 1e-2 * length;
            });
        std::transform(
            meshLengths.begin(),
            meshLengths.end(),
            rootLengths.begin(),
            [] (Scalar length) -> Scalar {
                return (1 + 2e-2) * length;
            });
        root = geo::AABBoxNode<Cell>(rootBase, rootLengths, nullptr);

        for (auto& rCell : cells) {
            root.insert(&rCell);
        }

        root.partition(targetLeafWidth, maxTreeDepth);
        //root.shrink();

        return BVH::flatten(
            root,
            [&cells] (Ref<const Cell> rCellData) -> unsigned {
                return std::distance<Ptr<const Cell>>(cells.data(), &rCellData);},
            std::allocator<std::byte>());
}


} // namespace cie::fem
