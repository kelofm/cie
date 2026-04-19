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

// --- STL Includes ---
#include <span>
#include <memory>


namespace cie::fem {


template <
    unsigned Dim,
    maths::Expression TIntegrand,
    class TQuadraturePointData = void>
class IntegrandProcessor {
public:
    static inline constexpr unsigned Dimension = Dim;

    struct Properties {
        std::optional<std::size_t>  integrandBatchSize;
        std::optional<std::uint8_t> verbosity;
    };

    IntegrandProcessor();

    virtual ~IntegrandProcessor();

    template <
        concepts::Iterator TCellIt,
        QuadratureRuleFactoryLike<
            typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>::Data,
            TQuadraturePointData
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
    void process(
        TCellIt itCellBegin,
        TCellIt itCellEnd,
        Ref<const TQuadratureRuleFactory> rQuadratureRuleFactory,
        TIntegrandFactory&& rIntegrandFactory,
        TIntegralSink&& rIntegralSink,
        Ref<const Properties> rExecutionProperties);

    template <
        concepts::Iterator TCellIt,
        QuadratureRuleFactoryLike<
            typename std::remove_const_t<typename std::iterator_traits<TCellIt>::value_type>,
            TQuadraturePointData
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
    void process(
        TCellIt itCellBegin,
        TCellIt itCellEnd,
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

private:
    template <
        class TCell,
        class TCellGetter,
        class TCellIt,
        class TQuadratureRuleFactory,
        class TIntegrandFactory,
        class TIntegralSink>
    void processImpl(
        TCellGetter&& rCellGetter,
        TCellIt itCellBegin,
        TCellIt itCellEnd,
        TQuadratureRuleFactory&& rQuadratureRuleFactory,
        TIntegrandFactory&& rIntegrandFactory,
        TIntegralSink&& rIntegralSink,
        Ref<const Properties> rExecutionProperties);
}; // class IntegrandProcessor


template <
    unsigned Dim,
    maths::StaticExpression TIntegrand,
    class TQuadraturePointData = void>
class ParallelIntegrandProcessor : public IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData> {
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
    maths::StaticExpression TIntegrand,
    class TQuadraturePointData = void>
class SYCLIntegrandProcessor : public IntegrandProcessor<Dim,TIntegrand,TQuadraturePointData> {
public:
    SYCLIntegrandProcessor(std::shared_ptr<sycl::queue> pQueue);

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
