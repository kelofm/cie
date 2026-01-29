#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/concurrency/inc/sycl.hpp"

// --- STL Includes ---
#include <span>
#include <memory>


namespace cie::fem {


template <
    unsigned Dim,
    maths::StaticExpression TIntegrand,
    class TQuadraturePointData = void>
class IntegrandProcessor {
public:
    static inline constexpr unsigned Dimension = Dim;

    struct Properties {
        std::optional<std::size_t> integrandBatchSize;
        std::optional<std::size_t> integrandsPerItem;
    };

    IntegrandProcessor();

    virtual ~IntegrandProcessor();

    template <
        GraphLike TMesh,
        QuadratureRuleFactoryLike<
            TMesh,
            typename TMesh::Vertex::Data,
            TQuadraturePointData
        > TQuadratureRuleFactory,
        concepts::FunctionWithSignature<
            TIntegrand,
            Ref<const TMesh>,
            Ref<const typename TMesh::Vertex::Data>
        > TIntegrandFactory,
        concepts::FunctionWithSignature<
            void,
            std::span<const VertexID>,
            std::span<const typename TIntegrand::Value>
        > TIntegralSink
    > void integrate(
        Ref<const TMesh> rMesh,
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
