#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/concurrency/inc/sycl.hpp"
#include "packages/io/inc/Traits.hpp"
#include "packages/stl_extension/inc/TrackedAllocator.hpp"

// --- STL Includes ---
#include <array>
#include <span>
#include <memory>
#include <optional>


namespace cie::fem {


template <
    unsigned Dim,
    maths::SliceableExpression TIntegrand,
    class TQuadraturePointData = void>
class IntegrandProcessor {
public:
    static inline constexpr unsigned Dimension = Dim;

    using ContextTiling = Tiling<
        TIntegrand::Context::Dimension,
        typename TIntegrand::Context::Value>;

    struct Properties {
        std::optional<std::size_t>  integrandBatchSize;
        std::optional<std::uint8_t> verbosity;
        /// Maximum tile extent for each Context axis; an empty optional means unbounded.
        std::array<std::optional<std::size_t>,TIntegrand::Context::Dimension> maxTileSizes {};
    };

    IntegrandProcessor();

    explicit IntegrandProcessor(RightRef<std::unique_ptr<ContextTiling>> rpTiling);

    virtual ~IntegrandProcessor();

    template <
        CellLike TCell,
        QuadratureRuleFactoryLike<
            TCell,
            TQuadraturePointData
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
    /// @details If a cell has an `activeContext()` method returning either Context or InContext,
    ///          that context is tiled. Otherwise the complete context bounds are used. Each input
    ///          context is tiled independently and all tile combinations are evaluated.
    void process(
        std::span<const TCell> cells,
        Ref<const TQuadratureRuleFactory> rQuadratureRuleFactory,
        TIntegrandFactory&& rIntegrandFactory,
        TIntegralSink&& rIntegralSink,
        Ref<const Properties> rExecutionProperties);

    virtual std::unique_ptr<Properties> makeDefaultProperties() const;

protected:
    virtual void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const Properties> rExecutionProperties);

    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class IntegrandProcessor


template <
    unsigned Dim,
    maths::StaticExpression TIntegrand,
    class TQuadraturePointData = void>
requires maths::SliceableExpression<TIntegrand>
class ParallelIntegrandProcessor
    : public IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData> {
public:
    ParallelIntegrandProcessor(Ref<mp::ThreadPoolBase> rThreads);

private:
    using Base = IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData>;

    void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const typename Base::Properties> rExecutionProperties) override;

    Ptr<mp::ThreadPoolBase> _pThreads;
}; // class ParallelIntegrandProcessor


#ifdef CIE_ENABLE_SYCL
template <
    unsigned Dim,
    maths::SliceableExpression TIntegrand,
    class TQuadraturePointData = void>
class SYCLIntegrandProcessor : public IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData> {
public:
    /// @brief Construct a processor evaluating integrands in SYCL kernels.
    /// @param pQueue Queue to submit the kernels to.
    /// @param pIntegrandMemoryStats Statistics object shared with the allocators
    ///        the integrand factory constructs the integrands with. It is used to
    ///        measure the dynamic memory footprint of each integrand, based on
    ///        which device memory is budgeted. Required for (and only used with)
    ///        dynamically sized integrands; the stats must not be shared with
    ///        integrand construction happening concurrently on other threads.
    SYCLIntegrandProcessor(
        std::shared_ptr<sycl::queue> pQueue,
        Ptr<AllocatorStats> pIntegrandMemoryStats = nullptr);

    ~SYCLIntegrandProcessor();

protected:
    using Base = IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData>;

    void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const typename Base::Properties> rExecutionProperties) override;

    struct Impl;
    std::unique_ptr<Impl> _pSYCLImpl;
}; // class SYCLIntegrandProcessor
#endif // CIE_ENABLE_SYCL


} // namespace cie::fem


#include "packages/utilities/impl/IntegrandProcessor_impl.hpp"
