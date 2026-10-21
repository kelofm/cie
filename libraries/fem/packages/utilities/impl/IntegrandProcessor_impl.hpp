#pragma once

// help the language server
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <format>
#include <limits>
#include <numeric>
#include <tuple>
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
    };

    std::vector<Cell> cells;
    std::vector<Tile> tiles;
    std::vector<typename OutContext::Value> outputMap;
    std::vector<QPoint> quadraturePoints;
    std::unique_ptr<typename IntegrandProcessor<Dim,TIntegrand,TQD>::ContextTiling> pTiling;

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
                const std::size_t cellQuadraturePointCount) -> void {
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
                        .iWorkEnd = 0ul});
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
                        appendCellWork(
                            rCell,
                            rIntegrandFactory(rCell),
                            iQuadraturePointBegin,
                            newPointCount);
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


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
struct SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::Impl {
    using Base = IntegrandProcessor<Dim,TIntegrand,TQD>;
    using Cell = typename Base::Impl::Cell;
    using Tile = typename Base::Impl::Tile;
    using QPoint = typename Base::Impl::QPoint;

    std::shared_ptr<sycl::queue> pQueue;

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

            std::size_t iCellBegin = 0ul;
            std::size_t iCellEnd = cellCount;
            while (iCellBegin < iCellEnd) {
                const std::size_t iCellMiddle = iCellBegin + (iCellEnd - iCellBegin) / 2ul;
                if (iWork < pCells[iCellMiddle].iWorkEnd)
                    iCellEnd = iCellMiddle;
                else
                    iCellBegin = iCellMiddle + 1ul;
            }
            if (cellCount <= iCellBegin) return;

            const Cell& rCell = pCells[iCellBegin];
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
};


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::SYCLIntegrandProcessor(
    std::shared_ptr<sycl::queue> pQueue)
    : _pSYCLImpl(new Impl) {
        _pSYCLImpl->pQueue = std::move(pQueue);
}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::~SYCLIntegrandProcessor() = default;


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
requires maths::SliceableExpression<TIntegrand>
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
                maxBufferSize = std::max(maxBufferSize, rCell.integrand.bufferSize(rTile.inContext));
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
        const sycl::nd_range<1> range(workItemsPerGroup * groupCount, workItemsPerGroup);
        const std::size_t bufferStride = maxBufferSize + maxResultSize;

        auto pDeviceCells = makeDeviceMemory<typename Impl::Cell>(rCells.size(), rQueue);
        auto pDeviceTiles = makeDeviceMemory<typename Impl::Tile>(rTiles.size(), rQueue);
        auto pDeviceOutputMap = makeDeviceMemory<typename TIntegrand::OutContext::Value>(rOutputMap.size(), rQueue);
        auto pDeviceQuadraturePoints = makeDeviceMemory<typename Impl::QPoint>(rQuadraturePoints.size(), rQueue);
        auto pDeviceOutput = makeDeviceMemory<typename TIntegrand::Value>(output.size(), rQueue);
        auto pDeviceBuffer = makeDeviceMemory<typename TIntegrand::Value>(workItemCount * bufferStride, rQueue);

        rQueue.copy(rCells.data(), pDeviceCells.get(), rCells.size());
        rQueue.copy(rTiles.data(), pDeviceTiles.get(), rTiles.size());
        rQueue.copy(rOutputMap.data(), pDeviceOutputMap.get(), rOutputMap.size());
        rQueue.copy(rQuadraturePoints.data(), pDeviceQuadraturePoints.get(), rQuadraturePoints.size());
        rQueue.fill(pDeviceOutput.get(), static_cast<typename TIntegrand::Value>(0), output.size());
        rQueue.wait_and_throw();

        rQueue.submit([&] (sycl::handler& rHandler) {
            rHandler.parallel_for<typename Impl::Kernel>(
                range,
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
        rQueue.copy(pDeviceOutput.get(), output.data(), output.size()).wait_and_throw();
}


#endif // CIE_ENABLE_SYCL


} // namespace cie::fem
