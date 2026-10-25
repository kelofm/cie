#pragma once

// help the language server
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/io/inc/Serializer.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/stl_extension/inc/LinearAllocator.hpp"
#include "packages/stl_extension/inc/TrackedAllocator.hpp"

// --- STL Includes ---
#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <format>
#include <limits>
#include <new>
#include <numeric>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <vector>


namespace cie::fem {


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
struct IntegrandProcessor<Dim,TIntegrand,TQD>::Impl {
    using QPoint = QuadraturePoint<Dimension,typename TIntegrand::Value,TQD>;
    using Context = typename TIntegrand::Context;
    using InContext = typename TIntegrand::InContext;
    using OutContext = typename TIntegrand::OutContext;

    static constexpr std::size_t InputCount = std::tuple_size_v<InContext>;

    struct Tile {
        InContext inContext;
        std::size_t iOutputMapBegin;
        std::size_t outputCount;
    };

    struct Cell {
        TIntegrand integrand;
        VertexID cellID;
        std::size_t iQuadraturePointBegin;
        std::size_t quadraturePointCount;
        std::size_t iTileBegin;
        std::size_t tileCount;
        std::size_t iOutputBegin;
        std::size_t outputSize;
        std::size_t iWorkBegin;
        std::size_t iWorkEnd;
        /// Dynamic memory allocated while constructing the integrand, measured
        /// through the integrand allocators' shared AllocatorStats (0 if untracked).
        std::size_t integrandMemory;
        /// Number of dynamic allocations made while constructing the integrand.
        std::size_t integrandAllocationCount;
    };

    std::vector<Cell> cells;
    std::vector<Tile> tiles;
    std::vector<typename OutContext::Value> outputMap;
    std::vector<QPoint> quadraturePoints;
    std::unique_ptr<typename IntegrandProcessor<Dim,TIntegrand,TQD>::ContextTiling> pTiling;
    /// Shared statistics object measuring the dynamic memory of integrands
    /// constructed by the integrand factory; null if measurement is disabled.
    Ptr<AllocatorStats> pIntegrandStats = nullptr;

    std::size_t assignWorkRanges() {
        std::size_t iWorkEnd = 0ul;
        for (auto& rCell : cells) {
            rCell.iWorkBegin = iWorkEnd;
            CIE_CHECK(
                rCell.tileCount == 0ul
                    || rCell.quadraturePointCount <= std::numeric_limits<std::size_t>::max() / rCell.tileCount,
                "Tiled integrand work count overflow")
            const std::size_t cellWorkCount = rCell.quadraturePointCount * rCell.tileCount;
            CIE_CHECK(
                cellWorkCount <= std::numeric_limits<std::size_t>::max() - iWorkEnd,
                "Tiled integrand work count overflow")
            iWorkEnd += cellWorkCount;
            rCell.iWorkEnd = iWorkEnd;
        }
        return iWorkEnd;
    }
}; // struct IntegrandProcessor::Impl


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
IntegrandProcessor<Dim,TIntegrand,TQD>::IntegrandProcessor()
    : IntegrandProcessor([] {
        using Context = typename TIntegrand::Context;
        std::array<typename Context::Value,Context::Dimension> tileShape;
        tileShape.fill(std::numeric_limits<typename Context::Value>::max());
        return std::make_unique<ContextTiling>(tileShape);
    }())
{}


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
IntegrandProcessor<Dim,TIntegrand,TQD>::IntegrandProcessor(
    RightRef<std::unique_ptr<ContextTiling>> rpTiling)
    : _pImpl(new Impl) {
        CIE_CHECK(rpTiling, "IntegrandProcessor requires a non-null tiling policy")
        _pImpl->pTiling = std::move(rpTiling);
}


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
IntegrandProcessor<Dim,TIntegrand,TQD>::~IntegrandProcessor() = default;


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
template <
    CellLike TCell,
    QuadratureRuleFactoryLike<
        TCell,
        TQD
    > TQuadratureRuleFactory,
    concepts::FunctionWithSignature<
        TIntegrand,
        Ref<const TCell>
    > TIntegrandFactory,
    concepts::FunctionWithSignature<
        void,
        std::span<const VertexID>,
        std::span<const typename TIntegrand::Value>
    > TIntegralSink
>
void IntegrandProcessor<Dim,TIntegrand,TQD>::process(
    std::span<const TCell> cells,
    Ref<const TQuadratureRuleFactory> rQuadratureRuleFactory,
    TIntegrandFactory&& rIntegrandFactory,
    TIntegralSink&& rIntegralSink,
    Ref<const Properties> rExecutionProperties) {
        using Context = typename TIntegrand::Context;
        using InContext = typename TIntegrand::InContext;
        using OutContext = typename TIntegrand::OutContext;

        CIE_BEGIN_EXCEPTION_TRACING
            auto pDefaultProperties = this->makeDefaultProperties();
            const std::size_t batchSize = rExecutionProperties.integrandBatchSize.has_value()
                ? rExecutionProperties.integrandBatchSize.value()
                : pDefaultProperties->integrandBatchSize.value();
            CIE_CHECK(0ul < batchSize, "IntegrandProcessor batch size must be positive")

            const ContextTiling* pTiling = _pImpl->pTiling.get();
            std::unique_ptr<ContextTiling> pConfiguredTiling;
            const bool hasConfiguredTileSize = std::any_of(
                rExecutionProperties.maxTileSizes.begin(),
                rExecutionProperties.maxTileSizes.end(),
                [] (const auto& rSize) { return rSize.has_value(); });
            if (hasConfiguredTileSize) {
                std::array<typename Context::Value,Context::Dimension> tileShape{};
                for (unsigned iDimension=0u; iDimension<Context::Dimension; ++iDimension) {
                    const std::size_t tileSize = rExecutionProperties.maxTileSizes[iDimension]
                        .value_or(static_cast<std::size_t>(std::numeric_limits<typename Context::Value>::max()));
                    CIE_CHECK(
                        0ul < tileSize
                            && tileSize <= static_cast<std::size_t>(std::numeric_limits<typename Context::Value>::max()),
                        "Configured tile sizes must fit the integrand context index type")
                    tileShape[iDimension] = static_cast<typename Context::Value>(tileSize);
                }
                pConfiguredTiling = std::make_unique<ContextTiling>(tileShape);
                pTiling = pConfiguredTiling.get();
            }

            auto& rQuadraturePoints = _pImpl->quadraturePoints;
            auto& rCells = _pImpl->cells;
            auto& rTiles = _pImpl->tiles;
            auto& rOutputMap = _pImpl->outputMap;
            rCells.clear();
            rTiles.clear();
            rOutputMap.clear();
            rQuadraturePoints.resize(batchSize);

            std::vector<VertexID> cellIDs;
            std::vector<typename TIntegrand::Value> output;
            std::size_t quadraturePointCount = 0ul;
            std::size_t outputSize = 0ul;

            const auto appendCellWork = [&] (
                Ref<const TCell> rCell,
                TIntegrand&& rIntegrand,
                const std::size_t iQuadraturePointBegin,
                const std::size_t cellQuadraturePointCount,
                const std::size_t integrandMemory,
                const std::size_t integrandAllocationCount) -> void {
                    // Create the complete input context used to determine the dense output shape.
                    const auto contextBounds = rIntegrand.contextBounds();
                    Context completeContext;
                    for (unsigned iDimension=0u; iDimension<Context::Dimension; ++iDimension) {
                        completeContext.data()[2 * iDimension] = 0;
                        completeContext.data()[2 * iDimension + 1] = contextBounds[iDimension];
                    }
                    InContext completeInContext;
                    completeInContext.fill(completeContext);
                    const OutContext completeOutContext = rIntegrand.outContext(completeInContext);
                    const std::size_t cellOutputSize = rIntegrand.size(completeInContext);
                    CIE_CHECK(
                        completeOutContext.count() == cellOutputSize,
                        "A complete integrand context must map to its complete output")

                    std::array<typename OutContext::Value,OutContext::Dimension> outputBounds{};
                    for (unsigned iDimension=0u; iDimension<OutContext::Dimension; ++iDimension) {
                        CIE_CHECK(
                            completeOutContext.data()[2 * iDimension] == 0,
                            "Complete output contexts must be zero-based")
                        outputBounds[iDimension] = completeOutContext.data()[2 * iDimension + 1];
                    }

                    // A cell may provide either one active slab (shared by each operand) or the
                    // full operand-context array. If it does not expose activeContext(), preserve
                    // the legacy behavior and activate the complete ansatz space.
                    InContext activeInContext = completeInContext;
                    if constexpr (requires { rCell.activeContext(); }) {
                        using ActiveContext = std::remove_cvref_t<decltype(rCell.activeContext())>;
                        if constexpr (std::is_same_v<ActiveContext,InContext>) {
                            activeInContext = rCell.activeContext();
                        } else if constexpr (std::is_same_v<ActiveContext,Context>) {
                            activeInContext.fill(rCell.activeContext());
                        } else {
                            static_assert(
                                std::is_same_v<ActiveContext,Context> || std::is_same_v<ActiveContext,InContext>,
                                "Cell::activeContext() must return the integrand Context or InContext");
                        }
                    }

                    // Tile each input independently and evaluate every combination. This is
                    // necessary for bilinear forms: a row tile paired with a different column
                    // tile contributes an off-diagonal block.
                    std::array<std::vector<Context>,Impl::InputCount> axisTiles;
                    for (std::size_t iInput=0ul; iInput<Impl::InputCount; ++iInput) {
                        const std::size_t tileCount = pTiling->count(activeInContext[iInput]);
                        axisTiles[iInput].resize(tileCount);
                        CIE_CHECK(
                            pTiling->tile(activeInContext[iInput], axisTiles[iInput]),
                            "Failed to tile an active input context")
                    }

                    const std::size_t iTileBegin = rTiles.size();
                    InContext currentInContext;
                    const auto appendTileCombination = [&] (auto&& rSelf, const std::size_t iInput) -> void {
                        if (iInput < Impl::InputCount) {
                            for (const Context& rContext : axisTiles[iInput]) {
                                currentInContext[iInput] = rContext;
                                rSelf(rSelf, iInput + 1ul);
                            }
                            return;
                        }

                        const OutContext outContext = rIntegrand.outContext(currentInContext);
                        const std::size_t tileOutputCount = rIntegrand.size(currentInContext);
                        CIE_CHECK(
                            outContext.count() == tileOutputCount,
                            "A sliced integrand result must match its output hyperslab")
                        const std::size_t iOutputMapBegin = rOutputMap.size();
                        rOutputMap.resize(iOutputMapBegin + tileOutputCount);
                        const auto localMap = std::span<typename OutContext::Value>(rOutputMap)
                            .subspan(iOutputMapBegin, tileOutputCount);
                        CIE_CHECK(
                            outContext.expand(outputBounds, localMap),
                            "A sliced output hyperslab lies outside the complete output")

                        rTiles.push_back(typename Impl::Tile {
                            .inContext = currentInContext,
                            .iOutputMapBegin = iOutputMapBegin,
                            .outputCount = tileOutputCount});
                    };
                    appendTileCombination(appendTileCombination, 0ul);

                    rCells.push_back(typename Impl::Cell {
                        .integrand = std::move(rIntegrand),
                        .cellID = rCell.id(),
                        .iQuadraturePointBegin = iQuadraturePointBegin,
                        .quadraturePointCount = cellQuadraturePointCount,
                        .iTileBegin = iTileBegin,
                        .tileCount = rTiles.size() - iTileBegin,
                        .iOutputBegin = outputSize,
                        .outputSize = cellOutputSize,
                        .iWorkBegin = 0ul,
                        .iWorkEnd = 0ul,
                        .integrandMemory = integrandMemory,
                        .integrandAllocationCount = integrandAllocationCount});
                    outputSize += cellOutputSize;
            };

            const auto flushBatch = [&] () -> void {
                if (quadraturePointCount == 0ul) return;

                rQuadraturePoints.resize(quadraturePointCount);
                output.resize(outputSize);
                this->execute(output, rExecutionProperties);

                cellIDs.resize(rCells.size());
                std::transform(
                    rCells.begin(),
                    rCells.end(),
                    cellIDs.begin(),
                    [] (Ref<const typename Impl::Cell> rCell) -> VertexID { return rCell.cellID; });
                rIntegralSink(cellIDs, output);

                rCells.clear();
                rTiles.clear();
                rOutputMap.clear();
                outputSize = 0ul;
                quadraturePointCount = 0ul;
                rQuadraturePoints.resize(batchSize);
            };

            for (std::size_t iCell=0ul; iCell<cells.size(); ++iCell) {
                Ref<const TCell> rCell = cells[iCell];
                auto quadraturePointFactory = rQuadratureRuleFactory(rCell);
                bool cellAddedToBatch = false;

                while (true) {
                    if (quadraturePointCount == batchSize) {
                        flushBatch();
                        cellAddedToBatch = false;
                    }

                    std::span<typename Impl::QPoint> availableQuadraturePoints(
                        rQuadraturePoints.data() + quadraturePointCount,
                        batchSize - quadraturePointCount);
                    const std::size_t newPointCount = quadraturePointFactory(availableQuadraturePoints);
                    CIE_CHECK(
                        newPointCount <= availableQuadraturePoints.size(),
                        "Quadrature rule factory returned more points than the available span")

                    if (newPointCount == 0ul) {
                        if (quadraturePointCount == batchSize) {
                            flushBatch();
                            cellAddedToBatch = false;
                            continue;
                        }
                        break;
                    }

                    if (!cellAddedToBatch) {
                        const std::size_t iQuadraturePointBegin = quadraturePointCount;
                        // Bracket construction with cumulative allocator
                        // statistics. Do not reset them: earlier integrands may
                        // still own allocations tracked by the same stats object.
                        const std::size_t allocatedBytesBefore = _pImpl->pIntegrandStats
                            ? _pImpl->pIntegrandStats->current() : 0ul;
                        const std::size_t allocationCountBefore = _pImpl->pIntegrandStats
                            ? _pImpl->pIntegrandStats->count() : 0ul;
                        TIntegrand constructedIntegrand = rIntegrandFactory(rCell);
                        const std::size_t allocatedBytesAfter = _pImpl->pIntegrandStats
                            ? _pImpl->pIntegrandStats->current() : 0ul;
                        const std::size_t integrandMemory = allocatedBytesAfter - allocatedBytesBefore;
                        const std::size_t integrandAllocationCount = _pImpl->pIntegrandStats
                            ? _pImpl->pIntegrandStats->count() - allocationCountBefore : 0ul;
                        appendCellWork(
                            rCell,
                            std::move(constructedIntegrand),
                            iQuadraturePointBegin,
                            newPointCount,
                            integrandMemory,
                            integrandAllocationCount);
                        cellAddedToBatch = true;
                    } else {
                        rCells.back().quadraturePointCount += newPointCount;
                    }
                    quadraturePointCount += newPointCount;
                }
            }
            flushBatch();
        CIE_END_EXCEPTION_TRACING
}


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
std::unique_ptr<typename IntegrandProcessor<Dim,TIntegrand,TQD>::Properties>
IntegrandProcessor<Dim,TIntegrand,TQD>::makeDefaultProperties() const {
    return std::make_unique<Properties>(Properties {
        .integrandBatchSize = 0x8000,
        .verbosity = 1,
        .maxTileSizes = {}
    });
}


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
void IntegrandProcessor<Dim,TIntegrand,TQD>::execute(
    std::span<typename TIntegrand::Value> output,
    Ref<const Properties> rExecutionProperties) {
        const auto timer = utils::LoggerSingleton::get().startTimer();
        const auto& rQuadraturePoints = _pImpl->quadraturePoints;
        CIE_CHECK(output.size() == std::accumulate(
            _pImpl->cells.begin(),
            _pImpl->cells.end(),
            std::size_t(0),
            [] (std::size_t sum, Ref<const typename Impl::Cell> rCell) { return sum + rCell.outputSize; }),
            "Output span does not match the cell work batch")
        std::fill(output.begin(), output.end(), static_cast<typename TIntegrand::Value>(0));

        std::size_t integrandCount = 0ul;
        std::vector<typename TIntegrand::Value> result;
        std::vector<typename TIntegrand::Value> buffer;
        for (const auto& rCell : _pImpl->cells) {
            integrandCount += rCell.tileCount;
            for (std::size_t iQuadraturePoint=rCell.iQuadraturePointBegin;
                 iQuadraturePoint<rCell.iQuadraturePointBegin + rCell.quadraturePointCount;
                 ++iQuadraturePoint) {
                for (std::size_t iTile=0ul; iTile<rCell.tileCount; ++iTile) {
                    const auto& rTile = _pImpl->tiles[rCell.iTileBegin + iTile];
                    result.resize(rTile.outputCount);
                    buffer.resize(rCell.integrand.bufferSize(rTile.inContext));
                    rQuadraturePoints[iQuadraturePoint].evaluate(
                        rCell.integrand,
                        result,
                        buffer,
                        rTile.inContext);
                    for (std::size_t iComponent=0ul; iComponent<rTile.outputCount; ++iComponent) {
                        output[rCell.iOutputBegin + _pImpl->outputMap[rTile.iOutputMapBegin + iComponent]] += result[iComponent];
                    }
                }
            }
        }

        if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
            utils::LoggerSingleton::get().logElapsed(
                std::format(
                    "evaluated {} tiled integrand(s) at {} quadrature point(s) in",
                    integrandCount,
                    rQuadraturePoints.size()),
                timer);
        }
}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
ParallelIntegrandProcessor<Dim,TIntegrand,TQD>::ParallelIntegrandProcessor(Ref<mp::ThreadPoolBase> rThreads)
    : _pThreads(&rThreads)
{}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
void ParallelIntegrandProcessor<Dim,TIntegrand,TQD>::execute(
    std::span<typename TIntegrand::Value> output,
    Ref<const typename Base::Properties> rExecutionProperties) {
        const auto timer = utils::LoggerSingleton::get().startTimer();
        const auto& rCells = this->_pImpl->cells;
        const auto& rTiles = this->_pImpl->tiles;
        const auto& rQuadraturePoints = this->_pImpl->quadraturePoints;
        const std::size_t workCount = this->_pImpl->assignWorkRanges();
        std::fill(output.begin(), output.end(), static_cast<typename TIntegrand::Value>(0));
        if (workCount == 0ul) return;

        mp::ParallelFor<std::size_t>(*_pThreads)
            .firstPrivate(std::vector<typename TIntegrand::Value>(), std::vector<typename TIntegrand::Value>())
            .execute(
                workCount,
                [
                    output,
                    pCells = rCells.data(),
                    pCellEnd = rCells.data() + rCells.size(),
                    pTiles = rTiles.data(),
                    pOutputMap = this->_pImpl->outputMap.data(),
                    pQuadraturePoints = rQuadraturePoints.data()
                ] (
                    std::size_t iWork,
                    Ref<std::vector<typename TIntegrand::Value>> rResult,
                    Ref<std::vector<typename TIntegrand::Value>> rBuffer) {
                        const auto pCell = std::upper_bound(
                            pCells,
                            pCellEnd,
                            iWork,
                            [] (std::size_t iWork, Ref<const typename Base::Impl::Cell> rCell) {
                                return iWork < rCell.iWorkEnd;
                            });
                        const auto& rCell = *pCell;
                        const std::size_t iLocalWork = iWork - rCell.iWorkBegin;
                        const std::size_t iTile = iLocalWork % rCell.tileCount;
                        const std::size_t iQuadraturePoint =
                            rCell.iQuadraturePointBegin + iLocalWork / rCell.tileCount;
                        const auto& rTile = pTiles[rCell.iTileBegin + iTile];
                        rResult.resize(rTile.outputCount);
                        rBuffer.resize(rCell.integrand.bufferSize(rTile.inContext));
                        pQuadraturePoints[iQuadraturePoint].evaluate(
                            rCell.integrand,
                            rResult,
                            rBuffer,
                            rTile.inContext);
                        for (std::size_t iComponent=0ul; iComponent<rTile.outputCount; ++iComponent) {
                            auto& rOutput = output[
                                rCell.iOutputBegin
                                + pOutputMap[rTile.iOutputMapBegin + iComponent]];
                            std::atomic_ref<typename TIntegrand::Value>(rOutput) += rResult[iComponent];
                        }
                });

        if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
            utils::LoggerSingleton::get().logElapsed(
                std::format(
                    "evaluated {} tiled integrand(s) at {} quadrature point(s) in",
                    workCount,
                    rQuadraturePoints.size()),
                timer);
        }
}


#ifdef CIE_ENABLE_SYCL


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
struct SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::Impl {
    using Base = IntegrandProcessor<Dim,TIntegrand,TQD>;
    using Cell = typename Base::Impl::Cell;
    using Tile = typename Base::Impl::Tile;
    using QPoint = typename Base::Impl::QPoint;
    using OutContext = typename TIntegrand::OutContext;

    /// @brief The integrand type with its dynamically sized components rebound to a
    ///        linear allocator, enabling device-side deserialization into a memory pool.
    using DeviceIntegrand = typename TIntegrand::template Rebind<LinearAllocator,std::byte>;

    std::shared_ptr<sycl::queue> pQueue;

    /// @brief Evaluation kernel for trivially copyable integrands, which are
    ///        embedded in the cell records moved to the device with a single copy.
    struct Kernel {
        Ptr<const Cell> pCells;
        Ptr<const Tile> pTiles;
        Ptr<const typename OutContext::Value> pOutputMap;
        Ptr<const QPoint> pQuadraturePoints;
        Ptr<typename TIntegrand::Value> pOutput;
        Ptr<typename TIntegrand::Value> pBuffer;
        std::size_t workItemCount;
        std::size_t cellCount;
        std::size_t maxBufferSize;
        std::size_t bufferStride;

        void operator()(sycl::nd_item<1> item) const noexcept {
            const std::size_t iWork = item.get_global_linear_id();
            if (workItemCount <= iWork) return;

            const auto pCell = std::upper_bound(
                pCells,
                pCells + cellCount,
                iWork,
                [] (std::size_t iWork, const Cell& rCell) { return iWork < rCell.iWorkEnd; });
            if (pCell == pCells + cellCount) return;

            const Cell& rCell = *pCell;
            const std::size_t iLocalWork = iWork - rCell.iWorkBegin;
            const std::size_t iTile = iLocalWork % rCell.tileCount;
            const std::size_t iQuadraturePoint =
                rCell.iQuadraturePointBegin + iLocalWork / rCell.tileCount;
            const Tile& rTile = pTiles[rCell.iTileBegin + iTile];
            const std::size_t iBufferBegin = iWork * bufferStride;
            const std::span<typename TIntegrand::Value> buffer(
                pBuffer + iBufferBegin,
                rCell.integrand.bufferSize(rTile.inContext));
            const std::span<typename TIntegrand::Value> result(
                pBuffer + iBufferBegin + maxBufferSize,
                rTile.outputCount);
            pQuadraturePoints[iQuadraturePoint].evaluate(
                rCell.integrand,
                result,
                buffer,
                rTile.inContext);

            for (std::size_t iComponent=0ul; iComponent<rTile.outputCount; ++iComponent) {
                auto& rOutput = pOutput[
                    rCell.iOutputBegin
                    + pOutputMap[rTile.iOutputMapBegin + iComponent]];
                sycl::atomic_ref<
                    typename TIntegrand::Value,
                    sycl::memory_order::relaxed,
                    sycl::memory_scope::device,
                    sycl::access::address_space::global_space
                >(rOutput) += result[iComponent];
            }
        }
    };

    /// @brief Cell record for dynamically sized integrands, which live at an
    ///        offset in the device-side memory pool.
    struct DynamicCell {
        std::size_t iQuadraturePointBegin;
        std::size_t quadraturePointCount;
        std::size_t iTileBegin;
        std::size_t tileCount;
        std::size_t iOutputBegin;
        std::size_t outputSize;
        std::size_t iWorkBegin;
        std::size_t iWorkEnd;
        /// Range of the integrand's serialized representation.
        std::size_t iSerializedBegin;
        std::size_t serializedSize;
        /// Range of the integrand's block in the device memory pool: the object
        /// storage is followed by persistent arena metadata and its allocation buffer.
        std::size_t iPoolObjectBegin;
        std::size_t iPoolArenaBegin;
        std::size_t iPoolBufferBegin;
        std::size_t arenaSize;
    }; // struct DynamicCell

    /// @brief Kernel deserializing one integrand per work item from its
    ///        serialized representation into the device memory pool.
    struct DeserializationKernel {
        Ptr<const DynamicCell> pCells;
        Ptr<std::byte> pPool;
        Ptr<const std::byte> pSerialized;
        std::size_t cellCount;

        void operator()(sycl::nd_item<1> item) const noexcept {
            const std::size_t iCell = item.get_global_linear_id();
            if (cellCount <= iCell) return;
            const DynamicCell& rCell = pCells[iCell];

            // Default-construct the integrand; its empty dynamic components hold
            // null-arena allocators that are replaced during deserialization.
            DeviceIntegrand* const pIntegrand
                = reinterpret_cast<DeviceIntegrand*>(pPool + rCell.iPoolObjectBegin);
            new (pIntegrand) DeviceIntegrand();

            // Keep the arena metadata in the pool as allocator instances retain
            // a pointer to it after this setup kernel completes.
            LinearArena* const pArena = reinterpret_cast<LinearArena*>(
                pPool + rCell.iPoolArenaBegin);
            new (pArena) LinearArena();
            pArena->buffer = std::span<std::byte>(
                pPool + rCell.iPoolBufferBegin,
                rCell.arenaSize);
            pArena->state = 0ul;
            LinearAllocator<std::byte> allocator(*pArena);

            cie::io::ByteDeserializer serialized(
                pSerialized + rCell.iSerializedBegin,
                rCell.serializedSize);
            cie::io::BinarySerializer::deserialize(
                serialized,
                *pIntegrand,
                allocator);
        }
    }; // struct DeserializationKernel

    /// @brief Evaluation kernel for dynamically sized integrands.
    struct DynamicKernel {
        Ptr<const DynamicCell> pCells;
        Ptr<const std::byte> pPool;
        Ptr<const Tile> pTiles;
        Ptr<const typename OutContext::Value> pOutputMap;
        Ptr<const QPoint> pQuadraturePoints;
        Ptr<typename TIntegrand::Value> pOutput;
        Ptr<typename TIntegrand::Value> pBuffer;
        std::size_t workItemCount;
        std::size_t cellCount;
        std::size_t maxBufferSize;
        std::size_t bufferStride;

        void operator()(sycl::nd_item<1> item) const noexcept {
            const std::size_t iWork = item.get_global_linear_id();
            if (workItemCount <= iWork) return;

            const auto pCell = std::upper_bound(
                pCells,
                pCells + cellCount,
                iWork,
                [] (std::size_t iWork, const DynamicCell& rCell) { return iWork < rCell.iWorkEnd; });
            if (pCell == pCells + cellCount) return;

            const DynamicCell& rCell = *pCell;
            const DeviceIntegrand& rIntegrand = *reinterpret_cast<const DeviceIntegrand*>(
                pPool + rCell.iPoolObjectBegin);
            const std::size_t iLocalWork = iWork - rCell.iWorkBegin;
            const std::size_t iTile = iLocalWork % rCell.tileCount;
            const std::size_t iQuadraturePoint =
                rCell.iQuadraturePointBegin + iLocalWork / rCell.tileCount;
            const Tile& rTile = pTiles[rCell.iTileBegin + iTile];
            const std::size_t iBufferBegin = iWork * bufferStride;
            const std::span<typename TIntegrand::Value> buffer(
                pBuffer + iBufferBegin,
                rIntegrand.bufferSize(rTile.inContext));
            const std::span<typename TIntegrand::Value> result(
                pBuffer + iBufferBegin + maxBufferSize,
                rTile.outputCount);
            pQuadraturePoints[iQuadraturePoint].evaluate(
                rIntegrand,
                result,
                buffer,
                rTile.inContext);

            for (std::size_t iComponent=0ul; iComponent<rTile.outputCount; ++iComponent) {
                auto& rOutput = pOutput[
                    rCell.iOutputBegin
                    + pOutputMap[rTile.iOutputMapBegin + iComponent]];
                sycl::atomic_ref<
                    typename TIntegrand::Value,
                    sycl::memory_order::relaxed,
                    sycl::memory_scope::device,
                    sycl::access::address_space::global_space
                >(rOutput) += result[iComponent];
            }
        }
    }; // struct DynamicKernel
};


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::SYCLIntegrandProcessor(
    std::shared_ptr<sycl::queue> pQueue,
    Ptr<AllocatorStats> pIntegrandMemoryStats)
    : _pSYCLImpl(new Impl) {
        _pSYCLImpl->pQueue = std::move(pQueue);
        this->_pImpl->pIntegrandStats = pIntegrandMemoryStats;
}


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::~SYCLIntegrandProcessor() = default;


template <unsigned Dim, maths::SliceableExpression TIntegrand, class TQD>
void SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::execute(
    std::span<typename TIntegrand::Value> output,
    Ref<const typename Base::Properties> rExecutionProperties) {
        Ref<sycl::queue> rQueue = *_pSYCLImpl->pQueue;
        const auto& rCells = this->_pImpl->cells;
        const auto& rTiles = this->_pImpl->tiles;
        const auto& rOutputMap = this->_pImpl->outputMap;
        const auto& rQuadraturePoints = this->_pImpl->quadraturePoints;
        const std::size_t workItemCount = this->_pImpl->assignWorkRanges();
        std::size_t maxBufferSize = 0ul;
        std::size_t maxResultSize = 0ul;
        for (const auto& rCell : rCells) {
            for (std::size_t iTile=0ul; iTile<rCell.tileCount; ++iTile) {
                const auto& rTile = rTiles[rCell.iTileBegin + iTile];
                maxBufferSize = std::max(maxBufferSize, static_cast<std::size_t>(rCell.integrand.bufferSize(rTile.inContext)));
                maxResultSize = std::max(maxResultSize, rTile.outputCount);
            }
        }
        std::fill(output.begin(), output.end(), static_cast<typename TIntegrand::Value>(0));
        if (workItemCount == 0ul) return;

        const std::string deviceName = rQueue.get_device().get_info<sycl::info::device::name>();
        std::string_view trimmedDeviceName = deviceName;
        trimmedDeviceName.remove_suffix(
            std::distance(
                deviceName.crbegin(),
                std::find_if(deviceName.crbegin(), deviceName.crend(), [](int c) { return !std::isspace(c); })));
        std::optional<utils::LogBlock> maybeLogBlock;
        if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value())
            maybeLogBlock.emplace(
                std::format(
                    "evaluate {} tiled work item(s) on {}",
                    workItemCount,
                    trimmedDeviceName),
                utils::LoggerSingleton::get());

        constexpr std::size_t kernelMaxWorkGroupSize = 256;
        const std::size_t workItemsPerGroup = std::min(workItemCount, kernelMaxWorkGroupSize);
        const std::size_t groupCount = (workItemCount + workItemsPerGroup - 1ul) / workItemsPerGroup;
        const sycl::nd_range<1> evaluationRange(workItemsPerGroup * groupCount, workItemsPerGroup);
        const std::size_t bufferStride = maxBufferSize + maxResultSize;

        auto pDeviceTiles = makeDeviceMemory<typename Impl::Tile>(rTiles.size(), rQueue);
        auto pDeviceOutputMap = makeDeviceMemory<typename TIntegrand::OutContext::Value>(rOutputMap.size(), rQueue);
        auto pDeviceQuadraturePoints = makeDeviceMemory<typename Impl::QPoint>(rQuadraturePoints.size(), rQueue);
        auto pDeviceOutput = makeDeviceMemory<typename TIntegrand::Value>(output.size(), rQueue);
        auto pDeviceBuffer = makeDeviceMemory<typename TIntegrand::Value>(workItemCount * bufferStride, rQueue);

        rQueue.copy(rTiles.data(), pDeviceTiles.get(), rTiles.size());
        rQueue.copy(rOutputMap.data(), pDeviceOutputMap.get(), rOutputMap.size());
        rQueue.copy(rQuadraturePoints.data(), pDeviceQuadraturePoints.get(), rQuadraturePoints.size());
        rQueue.fill(pDeviceOutput.get(), static_cast<typename TIntegrand::Value>(0), output.size());
        rQueue.wait_and_throw();

        if constexpr (std::is_trivially_copyable_v<TIntegrand>) {
            // Trivially copyable integrands are embedded in the cell records and
            // moved to the device with a single copy.
            auto pDeviceCells = makeDeviceMemory<typename Impl::Cell>(rCells.size(), rQueue);

            rQueue.copy(rCells.data(), pDeviceCells.get(), rCells.size());
            rQueue.wait_and_throw();

            rQueue.submit([&] (sycl::handler& rHandler) {
                rHandler.parallel_for<typename Impl::Kernel>(
                    evaluationRange,
                    typename Impl::Kernel {
                        .pCells = pDeviceCells.get(),
                        .pTiles = pDeviceTiles.get(),
                        .pOutputMap = pDeviceOutputMap.get(),
                        .pQuadraturePoints = pDeviceQuadraturePoints.get(),
                        .pOutput = pDeviceOutput.get(),
                        .pBuffer = pDeviceBuffer.get(),
                        .workItemCount = workItemCount,
                        .cellCount = rCells.size(),
                        .maxBufferSize = maxBufferSize,
                        .bufferStride = bufferStride});
            });
            rQueue.wait_and_throw();
        } else {
            // Dynamically sized integrands are serialized on the host, moved to
            // the device, then deserialized on the device into a memory pool from
            // which their dynamic components allocate through a linear allocator.
            using DeviceIntegrand = typename Impl::DeviceIntegrand;
            static_assert(
                BinarySerializable<TIntegrand>,
                "Dynamically sized integrands must be serializable");
            static_assert(
                std::is_default_constructible_v<DeviceIntegrand>,
                "Dynamically sized integrands must be default constructible "
                "after rebinding their allocators");
            if constexpr (!maths::StaticExpression<TIntegrand>) {
                CIE_CHECK(
                    this->_pImpl->pIntegrandStats != nullptr,
                    "SYCLIntegrandProcessor requires integrand allocator statistics "
                    "(AllocatorStats) for dynamically sized integrands")
            }

            // Serialize the integrands and lay out the device memory pool. Each
            // cell occupies one pool block: aligned object storage followed by its arena.
            constexpr std::size_t alignment = alignof(std::max_align_t);
            const auto alignUp = [] (const std::size_t size) {
                return (size + alignment - 1ul) / alignment * alignment;
            };
            std::vector<std::byte> serializedIntegrands;
            std::vector<typename Impl::DynamicCell> deviceCells(rCells.size());
            std::size_t poolSize = 0ul;
            for (std::size_t iCell=0ul; iCell<rCells.size(); ++iCell) {
                const auto& rCell = rCells[iCell];
                if constexpr (!maths::StaticExpression<TIntegrand>) {
                    CIE_CHECK(
                        this->_pImpl->pIntegrandStats != nullptr,
                        "Dynamic integrands must be constructed with tracked allocators "
                        "sharing the statistics object passed to SYCLIntegrandProcessor")
                }
                std::ostringstream stream;
                cie::io::BinarySerializer::serialize(stream, rCell.integrand);
                const std::string bytes = stream.str();

                auto& rDeviceCell = deviceCells[iCell];
                rDeviceCell.iSerializedBegin = serializedIntegrands.size();
                rDeviceCell.serializedSize = bytes.size();
                serializedIntegrands.insert(
                    serializedIntegrands.end(),
                    reinterpret_cast<Ptr<const std::byte>>(bytes.data()),
                    reinterpret_cast<Ptr<const std::byte>>(bytes.data()) + bytes.size());
                rDeviceCell.iPoolObjectBegin = poolSize;
                rDeviceCell.iPoolArenaBegin = poolSize + alignUp(sizeof(DeviceIntegrand));
                rDeviceCell.iPoolBufferBegin = rDeviceCell.iPoolArenaBegin + alignUp(sizeof(LinearArena));
                // The arena accommodates the integrand's tracked dynamic memory
                // plus the worst-case alignment padding of one allocation each.
                const std::size_t trackedArenaSize = rCell.integrandMemory
                    + (alignment - 1ul) * (rCell.integrandAllocationCount + 2ul);
                // Leave a small safety reserve for allocator/container behavior
                // that can differ between host construction and device-side
                // deserialization (for example, vector growth policies).
                rDeviceCell.arenaSize = std::max<std::size_t>(trackedArenaSize, 4096ul);
                poolSize = rDeviceCell.iPoolBufferBegin + alignUp(rDeviceCell.arenaSize);
                rDeviceCell.iQuadraturePointBegin = rCell.iQuadraturePointBegin;
                rDeviceCell.quadraturePointCount = rCell.quadraturePointCount;
                rDeviceCell.iTileBegin = rCell.iTileBegin;
                rDeviceCell.tileCount = rCell.tileCount;
                rDeviceCell.iOutputBegin = rCell.iOutputBegin;
                rDeviceCell.outputSize = rCell.outputSize;
                rDeviceCell.iWorkBegin = rCell.iWorkBegin;
                rDeviceCell.iWorkEnd = rCell.iWorkEnd;
            }

            auto pDeviceCells = makeDeviceMemory<typename Impl::DynamicCell>(deviceCells.size(), rQueue);
            auto pDevicePool = makeDeviceMemory<std::byte>(poolSize, rQueue);
            auto pDeviceSerialized = makeDeviceMemory<std::byte>(serializedIntegrands.size(), rQueue);

            rQueue.copy(deviceCells.data(), pDeviceCells.get(), deviceCells.size());
            rQueue.copy(serializedIntegrands.data(), pDeviceSerialized.get(), serializedIntegrands.size());
            rQueue.wait_and_throw();

            // Deserialize one integrand per device work item, with an isolated
            // arena block for each integrand.
            {
                const std::size_t deserializationWorkItemsPerGroup
                    = std::min(rCells.size(), kernelMaxWorkGroupSize);
                const std::size_t deserializationGroupCount
                    = (rCells.size() + deserializationWorkItemsPerGroup - 1ul)
                    / deserializationWorkItemsPerGroup;
                const sycl::nd_range<1> deserializationRange(
                    deserializationWorkItemsPerGroup * deserializationGroupCount,
                    deserializationWorkItemsPerGroup);
                rQueue.submit([&] (sycl::handler& rHandler) {
                    rHandler.parallel_for<typename Impl::DeserializationKernel>(
                        deserializationRange,
                        typename Impl::DeserializationKernel {
                            .pCells = pDeviceCells.get(),
                            .pPool = pDevicePool.get(),
                            .pSerialized = pDeviceSerialized.get(),
                            .cellCount = rCells.size()});
                });
                rQueue.wait_and_throw();
            }

            // Evaluate all dynamic integrands in parallel over the batch.
            rQueue.submit([&] (sycl::handler& rHandler) {
                rHandler.parallel_for<typename Impl::DynamicKernel>(
                    evaluationRange,
                    typename Impl::DynamicKernel {
                        .pCells = pDeviceCells.get(),
                        .pPool = pDevicePool.get(),
                        .pTiles = pDeviceTiles.get(),
                        .pOutputMap = pDeviceOutputMap.get(),
                        .pQuadraturePoints = pDeviceQuadraturePoints.get(),
                        .pOutput = pDeviceOutput.get(),
                        .pBuffer = pDeviceBuffer.get(),
                        .workItemCount = workItemCount,
                        .cellCount = rCells.size(),
                        .maxBufferSize = maxBufferSize,
                        .bufferStride = bufferStride});
            });
            rQueue.wait_and_throw();
        }
        rQueue.copy(pDeviceOutput.get(), output.data(), output.size()).wait_and_throw();
}


#endif // CIE_ENABLE_SYCL


} // namespace cie::fem
