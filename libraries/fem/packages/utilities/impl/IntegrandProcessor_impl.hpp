#pragma once

// help the language server
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"


namespace cie::fem {


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
struct IntegrandProcessor<TIntegrand,TQuadraturePointFactory>::Impl {
    using QPoint = QuadraturePoint<
        TQuadraturePointFactory::Dimension,
        typename TQuadraturePointFactory::Value>;

    class Extents {
    public:
        struct Value {
            VertexID cellID;
            std::size_t iQuadraturePointEnd;
            TIntegrand integrand;
        }; // struct Value

        Extents() {
            this->clear();
        }

        bool empty() const noexcept {
            return _data.size() == 1ul && _data.back().iQuadraturePointEnd == 0ul;
        }

        void push(VertexID cellID, RightRef<TIntegrand> rIntegrand) {
            _data.emplace_back(Value{
                .cellID = cellID,
                .iQuadraturePointEnd = _data.back().iQuadraturePoint,
                .integrand = std::move(rIntegrand)});
            std::swap(
                _data[_data.size() - 2],
                _data[_data.size() - 1]);
        }

        void extend(std::size_t count) noexcept {
            _data.back().iQuadraturePoint += count;
        }

        void clear() {
            _data.resize(1);
            _data.back().iQuadraturePoint = 0ul;
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

    OptionalRef<mp::ThreadPoolBase> maybeThreads;

    #ifdef CIE_ENABLE_SYCL
        DynamicArray<std::shared_ptr<sycl::queue>> queues;
    #endif

    Extents _extents;

    DynamicArray<QPoint> _quadraturePoints;
}; // struct IntegrandProcessor::Impl


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
IntegrandProcessor<TIntegrand,TQuadraturePointFactory>::IntegrandProcessor()
    : _pImpl(new Impl)
{}


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
IntegrandProcessor<TIntegrand,TQuadraturePointFactory>::IntegrandProcessor(OptionalRef<mp::ThreadPoolBase> rMaybeThreads)
    : IntegrandProcessor()
{
    _pImpl->maybeThreads = rMaybeThreads;
}


#ifdef CIE_ENABLE_SYCL
template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
IntegrandProcessor<TIntegrand,TQuadraturePointFactory>::IntegrandProcessor(
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads,
    std::span<const std::shared_ptr<sycl::queue>> queues)
    : IntegrandProcessor(rMaybeThreads)
{
    _pImpl->queues.insert(
        _pImpl->queues.end(),
        queues.begin(),
        queues.end());
}
#endif


template <
    maths::StaticExpression TIntegrand,
    QuadraturePointFactoryLike TQuadraturePointFactory>
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
> void IntegrandProcessor<TIntegrand,TQuadraturePointFactory>::integrate(
    Ref<const TMesh> rMesh,
    TIntegrandFactory&& rIntegrandFactory,
    TIntegralSink&& rIntegralSink)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_END_EXCEPTION_TRACING
}



} // namespace cie::fem
