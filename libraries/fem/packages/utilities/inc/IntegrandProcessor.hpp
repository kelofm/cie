#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/concurrency/inc/sycl.hpp"

// --- STL Includes ---
#include <span>
#include <memory>


namespace cie::fem {


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
class IntegrandProcessor {
    static_assert(
        std::is_same_v<
            typename TIntegrand::Value,
            typename TQuadraturePointFactory::Value>,
        "mismatched integrand and quadrature point value types");

public:
    struct Properties {
        std::optional<std::size_t> integrandBatchSize;
        std::optional<std::size_t> integrandsPerItem;
    };

    IntegrandProcessor();

    template <
        GraphLike TMesh,
        concepts::FunctionWithSignature<
            TIntegrand,
            Ref<const TMesh>,
            Ref<const typename TMesh::Vertex>
        > TIntegrandFactory,
        concepts::FunctionWithSignature<
            void,
            std::span<typename TIntegrand::Value>
        > TIntegralSink
    > void integrate(
        Ref<const TMesh> rMesh,
        TIntegrandFactory&& rIntegrandFactory,
        TIntegralSink&& rIntegralSink,
        Ref<const Properties> rExecutionProperties);

protected:
    virtual void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const Properties> rExecutionProperties);

    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class IntegrandProcessor


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
class ParallelIntegrandProcessor : public IntegrandProcessor<TIntegrand,TQuadraturePointFactory> {
public:
    ParallelIntegrandProcessor(Ref<mp::ThreadPoolBase> rThreads);

private:
    using Base = IntegrandProcessor<TIntegrand,TQuadraturePointFactory>;

    void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const typename Base::Properties> rExecutionProperties) override;

    Ptr<mp::ThreadPoolBase> _pThreads;
}; // class ParallelIntegrandProcessor


#ifdef CIE_ENABLE_SYCL
template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
class SYCLIntegrandProcessor : public IntegrandProcessor<TIntegrand,TQuadraturePointFactory> {
public:
    SYCLIntegrandProcessor(std::shared_ptr<sycl::queue> pQueue);

protected:
    using Base = IntegrandProcessor<TIntegrand,TQuadraturePointFactory>;

    void execute(
        std::span<typename TIntegrand::Value> output,
        Ref<const typename Base::Properties> rExecutionProperties) override;

    std::shared_ptr<sycl::queue> _pQueue;
}; // class SYCLIntegrandProcessor
#endif


} // namespace cie::fem
