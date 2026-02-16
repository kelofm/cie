#pragma once

// help the language server
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <array>
#include <algorithm>
#include <limits>
#include <format>


namespace cie::fem {


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
struct IntegrandProcessor<Dim,TIntegrand,TQD>::Impl {
    using QPoint = QuadraturePoint<Dimension,typename TIntegrand::Value,TQD>;

    struct StaticIntegrandSize {
        static constexpr std::size_t value() noexcept
        requires maths::StaticExpression<TIntegrand> {
            return TIntegrand::size();}

        static constexpr std::size_t value() noexcept
        requires (!maths::StaticExpression<TIntegrand>){
            return 0ul;}
    }; // struct StaticIntegrandSize

    class Extents {
    public:
        struct Value {
            VertexID cellID;
            std::size_t iQuadraturePointBegin;
            TIntegrand integrand;
        }; // struct Value

        Extents() {
            this->clear();
        }

        bool empty() const noexcept {
            return _data.size() == 1ul && _data.back().iQuadraturePointBegin == 0ul;
        }

        void push(VertexID cellID, RightRef<TIntegrand> rIntegrand) {
            if (1 < _data.size() && _data[_data.size() - 2].iQuadraturePointBegin == _data.back().iQuadraturePointBegin) [[unlikely]] {
                _data[_data.size() - 2] = Value {
                    .cellID = cellID,
                    .iQuadraturePointBegin = _data.back().iQuadraturePointBegin,
                    .integrand = std::move(rIntegrand)};
            } else {
                _data.emplace_back(Value{
                    .cellID = cellID,
                    .iQuadraturePointBegin = _data.back().iQuadraturePointBegin,
                    .integrand = std::move(rIntegrand)});
                std::swap(
                    _data[_data.size() - 2],
                    _data[_data.size() - 1]);
            }
        }

        void extend(std::size_t count) noexcept {
            _data.back().iQuadraturePointBegin += count;
        }

        void clear() {
            _data.resize(1);
            _data.back().iQuadraturePointBegin = 0ul;
            _data.back().cellID = std::numeric_limits<VertexID::Value>::max();
        }

        std::span<Value> get() noexcept {
            return _data;
        }

        std::span<const Value> get() const noexcept {
            return _data;
        }

    private:
        DynamicArray<Value> _data;
    }; // class Extents

    Extents extents;

    DynamicArray<QPoint> quadraturePoints;
}; // struct IntegrandProcessor::Impl


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
IntegrandProcessor<Dim,TIntegrand,TQD>::IntegrandProcessor()
    : _pImpl(new Impl)
{}


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
IntegrandProcessor<Dim,TIntegrand,TQD>::~IntegrandProcessor() = default;


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
template <
    concepts::Iterator TCellIt,
    QuadratureRuleFactoryLike<
        typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>::Data,
        TQD
    > TQuadratureRuleFactory,
    concepts::FunctionWithSignature<
        TIntegrand,
        Ref<const typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>::Data>
    > TIntegrandFactory,
    concepts::FunctionWithSignature<
        void,
        std::span<const VertexID>,
        std::span<const typename TIntegrand::Value>
    > TIntegralSink
>
requires CellLike<std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type::Data>>
void  IntegrandProcessor<Dim,TIntegrand,TQD>::process(
    TCellIt itCellBegin,
    TCellIt itCellEnd,
    Ref<const TQuadratureRuleFactory> rQuadratureRuleFactory,
    TIntegrandFactory&& rIntegrandFactory,
    TIntegralSink&& rIntegralSink,
    Ref<const Properties> rExecutionProperties)
{
    using TItValueType = typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>;
    this->processImpl<typename TItValueType::Data>(
        [] (const TItValueType& rItem) -> const typename TItValueType::Data& {return rItem.data();},
        itCellBegin,
        itCellEnd,
        rQuadratureRuleFactory,
        rIntegrandFactory,
        rIntegralSink,
        rExecutionProperties);
}


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
template <
    concepts::Iterator TCellIt,
    QuadratureRuleFactoryLike<
        typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>,
        TQD
    > TQuadratureRuleFactory,
    concepts::FunctionWithSignature<
        TIntegrand,
        Ref<const typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>>
    > TIntegrandFactory,
    concepts::FunctionWithSignature<
        void,
        std::span<const VertexID>,
        std::span<const typename TIntegrand::Value>
    > TIntegralSink
>
requires CellLike<std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>>
void  IntegrandProcessor<Dim,TIntegrand,TQD>::process(
    TCellIt itCellBegin,
    TCellIt itCellEnd,
    Ref<const TQuadratureRuleFactory> rQuadratureRuleFactory,
    TIntegrandFactory&& rIntegrandFactory,
    TIntegralSink&& rIntegralSink,
    Ref<const Properties> rExecutionProperties)
{
    using TItValueType = typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>;
    this->processImpl<TItValueType>(
        [] (const TItValueType& rItem) -> const TItValueType& {return rItem;},
        itCellBegin,
        itCellEnd,
        rQuadratureRuleFactory,
        rIntegrandFactory,
        rIntegralSink,
        rExecutionProperties);
}

template <unsigned Dim, maths::Expression TIntegrand, class TQD>
template <
    class TCell,
    class TCellGetter,
    class TCellIt,
    class TQuadratureRuleFactory,
    class TIntegrandFactory,
    class TIntegralSink>
void IntegrandProcessor<Dim,TIntegrand,TQD>::processImpl(
    TCellGetter&& rCellGetter,
    TCellIt itCellBegin,
    TCellIt itCellEnd,
    TQuadratureRuleFactory&& rQuadratureRuleFactory,
    TIntegrandFactory&& rIntegrandFactory,
    TIntegralSink&& rIntegralSink,
    Ref<const Properties> rExecutionProperties)
{
    CIE_BEGIN_EXCEPTION_TRACING
    // Parse execution properties.
    auto pDefaultProperties = this->makeDefaultProperties();
    const std::size_t batchSize = rExecutionProperties.integrandBatchSize.has_value()
        ? rExecutionProperties.integrandBatchSize.value()
        : pDefaultProperties->integrandBatchSize.value();

    // Allocate memory.
    _pImpl->quadraturePoints.resize(batchSize);
    DynamicArray<VertexID> cellIDs;
    DynamicArray<typename TIntegrand::Value> output;
    using TQuadraturePointFactory = std::invoke_result_t<
        TQuadratureRuleFactory,
        Ref<const TCell>>;
    TQuadraturePointFactory quadraturePointFactory;

    // Producer-consumer loop.
    // 0) Generate quadrature points.
    // 1) Evaluate integrands at quadrature points.
    // 2) Call the consumer with the cell-wise reduced results.
    for (auto itCell=itCellBegin; itCell!=itCellEnd; ++itCell) {
        Ref<const TCell> rCell = rCellGetter(*itCell);

        // Mark the begin of recording quadrature points for a new cell
        // and thus new integrand.
        _pImpl->extents.push(
            rCell.id(),
            rIntegrandFactory(rCell));

        // Construct a new quadrature point rule for the current cell.
        CIE_BEGIN_EXCEPTION_TRACING
        quadraturePointFactory = rQuadratureRuleFactory(rCell);
        CIE_END_EXCEPTION_TRACING

        // Generate quadrature points for the current cell,
        // filling the internal container.
        // If the internal container reaches capacity before
        // all quadrature points could be generated for the
        // current cell, evaluate  and clear the internal
        // structures, then continue generating quadrature
        // points.
        while (true) {
            // If an evaluation was issued and the internal structures cleared,
            // the current cell and its integrand must be pushed to the internal
            // structures.
            if (_pImpl->extents.empty()) {
                _pImpl->extents.push(
                    rCell.id(),
                    rIntegrandFactory(rCell));
            } // extents.empty()

            // Generate quadrature points into the remaining range
            // of available slots.
            const std::size_t iQuadraturePointBegin = _pImpl->extents.get().back().iQuadraturePointBegin;
            std::span<typename Impl::QPoint> availableQuadraturePointSpan(
                _pImpl->quadraturePoints.data() + iQuadraturePointBegin,
                _pImpl->quadraturePoints.data() + _pImpl->quadraturePoints.size());
            const std::size_t newQuadraturePointCount = quadraturePointFactory(availableQuadraturePointSpan);

            // Update the internal structures to account for the newly
            // generated quadrature points.
            _pImpl->extents.extend(newQuadraturePointCount);

            // Two reasons can cause no quadrature points to be generated.
            // 0) All quadrature points were generated for the current cell,
            //    in which case we can move on to the next cell without having
            //    to issue an evaluation.
            // 1) The internal quadrature point storage is full, in which case
            //    an evaluation must be issued and the internal structures cleared.
            //    Quadrature point generation can continue afterwards.
            if (newQuadraturePointCount == 0ul) {
                if (batchSize <= _pImpl->extents.get().back().iQuadraturePointBegin) {
                    // 1) The internal structures are full.
                    // => issue an evaluation.
                    const auto extentView = _pImpl->extents.get();
                    const std::size_t cellCount = extentView.size() - 1;
                    output.resize(cellCount * extentView.back().integrand.size());

                    CIE_BEGIN_EXCEPTION_TRACING
                        this->execute(
                            output,
                            rExecutionProperties);
                    CIE_END_EXCEPTION_TRACING

                    // Call the consumer with the reduced results.
                    cellIDs.resize(cellCount);
                    std::transform(
                        extentView.begin(),
                        extentView.begin() + cellCount,
                        cellIDs.begin(),
                        [](Ref<const typename Impl::Extents::Value> rEntry) -> VertexID {
                            return rEntry.cellID;
                        });
                    CIE_BEGIN_EXCEPTION_TRACING
                    rIntegralSink(
                        cellIDs,
                        output);
                    CIE_END_EXCEPTION_TRACING

                    // Clear the internal structures.
                    _pImpl->extents.clear();
                } /*if (internal structures are full)*/ else {
                    // 0) Finished generating quadrature points for the current cell.
                    break;
                } // if (internal structures are full) else
            } // if newQuadraturePointCount == 0ul
        } // while true
    } // for rCell in rMesh.vertices()

    // Issue an evaluation if there are any leftover quadrature points.
    if (!_pImpl->extents.empty() && _pImpl->extents.get().back().iQuadraturePointBegin != 0ul) {
        const auto extentView = _pImpl->extents.get();
        _pImpl->quadraturePoints.resize(extentView.back().iQuadraturePointBegin);
        const std::size_t cellCount = extentView.size() - 1;
        output.resize(cellCount * extentView.back().integrand.size());

        CIE_BEGIN_EXCEPTION_TRACING
            this->execute(
                output,
                rExecutionProperties);
        CIE_END_EXCEPTION_TRACING

        // Call the consumer with the reduced results.
        cellIDs.resize(cellCount);
        std::transform(
            extentView.begin(),
            extentView.end() - 1,
            cellIDs.begin(),
            [](Ref<const typename Impl::Extents::Value> rEntry) -> VertexID {
                return rEntry.cellID;
            });
        CIE_BEGIN_EXCEPTION_TRACING
        rIntegralSink(
            cellIDs,
            output);
        CIE_END_EXCEPTION_TRACING

        // Clear the internal structures.
        _pImpl->extents.clear();
    } // if !extents.empty()

    CIE_END_EXCEPTION_TRACING
}


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
std::unique_ptr<typename IntegrandProcessor<Dim,TIntegrand,TQD>::Properties>
IntegrandProcessor<Dim,TIntegrand,TQD>::makeDefaultProperties() const {
    return std::make_unique<Properties>(Properties {
        .integrandBatchSize = 0x8000,
        .integrandsPerItem = 0x8,
        .verbosity = 1
    });
}


template <unsigned Dim, maths::Expression TIntegrand, class TQD>
void IntegrandProcessor<Dim,TIntegrand,TQD>::execute(std::span<typename TIntegrand::Value> output,
                                                     Ref<const Properties> rExecutionProperties) {
    const auto timer = utils::LoggerSingleton::get().startTimer();

    const auto& rQuadraturePoints = _pImpl->quadraturePoints;
    const auto extentView = _pImpl->extents.get();
    if (rQuadraturePoints.empty()) return;

    std::conditional_t<
        maths::StaticExpression<TIntegrand>,
        std::array<typename TIntegrand::Value,Impl::StaticIntegrandSize::value()>,
        std::vector<typename TIntegrand::Value>
    > result;
    //std::array<typename TIntegrand::Value,TIntegrand::size()> result;

    if constexpr (!maths::StaticExpression<TIntegrand>) {
        result.resize(extentView.back().integrand.size());
    }

    for (std::size_t iExtent=0ul; iExtent<extentView.size()-1; ++iExtent) {
        Ref<const TIntegrand> rIntegrand = extentView[iExtent].integrand;
        const std::span<typename TIntegrand::Value> reducedOutput(
            output.data() + iExtent * result.size(),
            result.size());
        std::fill_n(
            reducedOutput.data(),
            result.size(),
            static_cast<typename TIntegrand::Value>(0));

        const std::size_t iQuadraturePointBegin = extentView[iExtent].iQuadraturePointBegin;
        const std::size_t iQuadraturePointEnd = extentView[iExtent + 1].iQuadraturePointBegin;
        for (std::size_t iQuadraturePoint=iQuadraturePointBegin; iQuadraturePoint<iQuadraturePointEnd; ++iQuadraturePoint) {
            rQuadraturePoints[iQuadraturePoint].evaluate(
                rIntegrand,
                result);
            const auto pLeftEnd = reducedOutput.data() + reducedOutput.size();
            auto pRight=result.data();
            for (auto pLeft=reducedOutput.data(); pLeft<pLeftEnd; ++pLeft, ++pRight) {
                *pLeft += *pRight;
            }
        } // for iQuadraturePoint in range(iQuadraturePointBegin, iQuadraturePointEnd)
    } // for iExtent in range(extentView.size())

    if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
        utils::LoggerSingleton::get().logElapsed(
            std::format(
                "evaluated {} integrand(s) at {} quadrature point(s) in",
                extentView.size() - 1,
                extentView.back().iQuadraturePointBegin),
            timer);
    }
}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
ParallelIntegrandProcessor<Dim,TIntegrand,TQD>::ParallelIntegrandProcessor(Ref<mp::ThreadPoolBase> rThreads)
    : _pThreads(&rThreads)
{}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
void ParallelIntegrandProcessor<Dim,TIntegrand,TQD>::execute(std::span<typename TIntegrand::Value> output,
                                                             Ref<const typename Base::Properties> rExecutionProperties) {
    const auto timer = utils::LoggerSingleton::get().startTimer();

    std::fill_n(
        output.data(),
        output.size(),
        static_cast<typename TIntegrand::Value>(0));

    const auto extentView = this->_pImpl->extents.get();
    mp::ParallelFor<std::size_t>(*_pThreads).operator()(
        this->_pImpl->quadraturePoints.size(),
        [
            output,
            pExtentBegin            = extentView.data(),
            pExtentEnd              = extentView.data() + extentView.size(),
            pQuadraturePointBegin   = this->_pImpl->quadraturePoints.data()
        ](std::size_t iQuadraturePoint){
            // Find which integrand the given quadrature point belongs to.
            Ptr<const typename Base::Impl::Extents::Value> pExtent = std::upper_bound(
                pExtentBegin,
                pExtentEnd,
                iQuadraturePoint,
                [](std::size_t iQuadraturePoint, Ref<const typename Base::Impl::Extents::Value> rExtent) -> bool {
                    return iQuadraturePoint < rExtent.iQuadraturePointBegin;
                });
            if (iQuadraturePoint < pExtent->iQuadraturePointBegin) --pExtent;

            // Evaluate the integrand.
            Ref<const typename Base::Impl::QPoint> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint];
            const TIntegrand integrand = pExtent->integrand;
            std::array<typename TIntegrand::Value,TIntegrand::size()> result;
            rQuadraturePoint.evaluate(integrand, result);

            // Reduce results.
            Ptr<typename TIntegrand::Value> pOutput = output.data() + std::distance<decltype(pExtent)>(pExtentBegin, pExtent) * TIntegrand::size();
            Ptr<const typename TIntegrand::Value> pOutputEnd = pOutput + TIntegrand::size();
            Ptr<const typename TIntegrand::Value> pResult = result.data();
            for (; pOutput<pOutputEnd; ++pOutput, ++pResult) {
                std::atomic_ref<typename TIntegrand::Value>(*pOutput) += *pResult;
            }
        });

    if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
        utils::LoggerSingleton::get().logElapsed(
            std::format(
                "evaluated {} integrand(s) at {} quadrature point(s) in",
                extentView.size() - 1,
                extentView.back().iQuadraturePointBegin),
            timer);
    }
}


#ifdef CIE_ENABLE_SYCL


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
struct SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::Impl {
    std::shared_ptr<sycl::queue> pQueue;
}; // struct SYCLIntegrandProcessor::Impl


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::SYCLIntegrandProcessor(
    std::shared_ptr<sycl::queue> pQueue)
    : _pSYCLImpl(new Impl)
{
    _pSYCLImpl->pQueue = pQueue;
}


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::~SYCLIntegrandProcessor() = default;


template <unsigned Dim, maths::StaticExpression TIntegrand, class TQD>
void SYCLIntegrandProcessor<Dim,TIntegrand,TQD>::execute(std::span<typename TIntegrand::Value> output,
                                                         Ref<const typename Base::Properties> rExecutionProperties) {
    const auto timer = utils::LoggerSingleton::get().startTimer();
    Ref<sycl::queue> rQueue = *_pSYCLImpl->pQueue;
    const auto extentView = this->_pImpl->extents.get();
    const auto& rQuadraturePoints = this->_pImpl->quadraturePoints;

    const std::string deviceName = rQueue.get_device().get_info<sycl::info::device::name>();
    std::string_view trimmedDeviceName = deviceName;
    trimmedDeviceName.remove_suffix(
        std::distance(
            deviceName.crbegin(),
            std::find_if(
                deviceName.crbegin(),
                deviceName.crend(),
                [](int c) {return !std::isspace(c);})));

    // Parse execution properties.
    auto pDefaultExecutionProperties = this->makeDefaultProperties();
    const std::size_t integrandsPerItem = rExecutionProperties.integrandsPerItem.has_value()
        ? rExecutionProperties.integrandsPerItem.value()
        : pDefaultExecutionProperties->integrandsPerItem.value();

    // Allocate memory on the device.
    if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
        utils::LoggerSingleton::get().log(std::format(
            "allocate {} bytes on {}",
            extentView.size() * sizeof(typename Base::Impl::Extents::Value)
                + rQuadraturePoints.size() * sizeof(typename Base::Impl::QPoint)
                + output.size() * sizeof(typename TIntegrand::Value),
            trimmedDeviceName));
    }
    auto pDeviceExtents = makeDeviceMemory<typename Base::Impl::Extents::Value>(
        extentView.size(),
        rQueue);
    auto pDeviceQuadraturePoints = makeDeviceMemory<typename Base::Impl::QPoint>(
        rQuadraturePoints.size(),
        rQueue);
    auto pDeviceOutput = makeDeviceMemory<typename TIntegrand::Value>(
        output.size(),
        rQueue);

    // Copy input data to the device.
    [[maybe_unused]] sycl::event extentCopyEvent = rQueue.copy(
        extentView.data(),
        pDeviceExtents.get(),
        extentView.size());
    [[maybe_unused]] sycl::event quadraturePointCopyEvent = rQueue.copy(
        rQuadraturePoints.data(),
        pDeviceQuadraturePoints.get(),
        rQuadraturePoints.size());

    // Initialize output values on the device side.
    [[maybe_unused]] sycl::event outputInitEvent = rQueue.fill(
        pDeviceOutput.get(),
        static_cast<typename TIntegrand::Value>(0),
        output.size());

    // Perform integration on the device.
    CIE_BEGIN_EXCEPTION_TRACING
    const std::size_t workItemCount = rQuadraturePoints.size() / integrandsPerItem + bool(rQuadraturePoints.size() % integrandsPerItem);
    const std::size_t itemsPerWorkGroup = std::min<std::size_t>(
        workItemCount,
        rQueue.get_device().get_info<sycl::info::device::max_work_group_size>());
    const std::size_t workGroupCount = workItemCount / itemsPerWorkGroup + bool(workItemCount % itemsPerWorkGroup);
    const auto range = sycl::nd_range<1>(
        itemsPerWorkGroup * workGroupCount,
        itemsPerWorkGroup);

    rQueue.wait_and_throw();
    rQueue.submit([&] (sycl::handler& rHandler) {
        rHandler.parallel_for(
            range,
            [
                pExtentBegin            = pDeviceExtents.get(),
                pExtentEnd              = pDeviceExtents.get() + extentView.size(),
                pQuadraturePointBegin   = pDeviceQuadraturePoints.get(),
                pOutputBegin            = pDeviceOutput.get(),
                quadraturePointCount    = rQuadraturePoints.size(),
                integrandsPerItem
            ] (sycl::nd_item<1> item) {
                // Get work item index.
                //const std::size_t iItem = item.get(0);
                const std::size_t iItem = item.get_global_linear_id();

                // Map to quadrature point range.
                const std::size_t iQuadraturePointBegin = std::min<std::size_t>(
                    iItem * integrandsPerItem,
                    quadraturePointCount);
                const std::size_t iQuadraturePointEnd   = std::min<std::size_t>(
                    iQuadraturePointBegin + integrandsPerItem,
                    quadraturePointCount);

                // Find which integrand the first quadrature point belongs to.
                Ptr<const typename Base::Impl::Extents::Value> pExtent = std::upper_bound(
                    pExtentBegin,
                    pExtentEnd,
                    iQuadraturePointBegin,
                    [](std::size_t iQuadraturePoint, Ref<const typename Base::Impl::Extents::Value> rExtent) -> bool {
                        return iQuadraturePoint < rExtent.iQuadraturePointBegin;
                    }) - 1;
                TIntegrand integrand = pExtent->integrand;
                Ptr<const typename Base::Impl::Extents::Value> pLastExtent = pExtent;
                std::array<typename TIntegrand::Value,TIntegrand::size()> results;

                // Loop through assigned quadrature points and evaluate their corresponding integrands.
                for (std::size_t iQuadraturePoint=iQuadraturePointBegin; iQuadraturePoint<iQuadraturePointEnd; ++iQuadraturePoint) {
                    // Find which integrand the current quadrature point belongs to.
                    for (; (pExtent + 1)->iQuadraturePointBegin <= iQuadraturePoint; ++pExtent) {}
                    if (pExtent != pLastExtent) integrand = pExtent->integrand;
                    pLastExtent = pExtent;

                    // Evaluate the integrand at the current quadrature point.
                    pQuadraturePointBegin[iQuadraturePoint].evaluate(
                        integrand,
                        results);

                    // Reduce the results.
                    std::span<typename TIntegrand::Value> output(
                        pOutputBegin + std::distance<decltype(pExtent)>(pExtentBegin, pExtent) * TIntegrand::size(),
                        TIntegrand::size());
                    for (std::size_t iComponent=0ul; iComponent<TIntegrand::size(); ++iComponent) {
                        sycl::atomic_ref<
                            typename TIntegrand::Value,
                            sycl::memory_order::relaxed,
                            sycl::memory_scope::device,
                            sycl::access::address_space::local_space
                        >(output[iComponent]) += results[iComponent];
                    } // for iComponent in range(TIntegrand::size())
                } // for iQuadraturePoint in range(iQuadraturePointBegin, iQuadraturePointEnd)
            });
        });
    rQueue.wait_and_throw();
    CIE_END_EXCEPTION_TRACING

    // Fetch results from the device.
    CIE_BEGIN_EXCEPTION_TRACING
    rQueue.copy(
        pDeviceOutput.get(),
        output.data(),
        output.size()).wait_and_throw();
    CIE_END_EXCEPTION_TRACING

    if (rExecutionProperties.verbosity && 3 <= rExecutionProperties.verbosity.value()) {
        utils::LoggerSingleton::get().logElapsed(
            std::format(
                "evaluated {} integrand(s) at {} quadrature point(s) on {} in",
                extentView.size() - 1,
                extentView.back().iQuadraturePointBegin,
                trimmedDeviceName),
            timer);
    }
}


#endif // CIE_ENABLE_SYCL



} // namespace cie::fem
