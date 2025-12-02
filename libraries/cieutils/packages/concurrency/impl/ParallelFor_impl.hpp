#ifndef CIE_UTILS_PARALLEL_FOR_IMPL_HPP
#define CIE_UTILS_PARALLEL_FOR_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/IndexPartitionFactory.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::mp {


template <concepts::Integer TIndex, class TStorage>
template <class ...TArgs>
requires (std::is_same_v<
            std::tuple<typename std::remove_reference_t<std::decay_t<TArgs>>...>,
            typename TStorage::StorageType // <== tuple of stored types
          >)
ParallelFor<TIndex,TStorage>::ParallelFor(Ref<ThreadPoolBase> r_pool, TArgs&&... rArgs) :
    _pool(r_pool, TStorage(std::forward<TArgs>(rArgs)...))
{
}


template <concepts::Integer TIndex, class TStorage>
template <class ...TArgs>
ParallelFor<TIndex,ThreadStorage<typename std::remove_reference_t<TArgs>...>>
ParallelFor<TIndex,TStorage>::firstPrivate(TArgs&&... rArgs)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return ParallelFor<
        TIndex,
        ThreadStorage<typename std::remove_reference_t<TArgs>...>
    >(_pool, std::forward<TArgs>(rArgs)...);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(Ref<const IndexPartitionFactory> rPartitionFactory,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(rPartitionFactory,
                  rFunction,
                  true);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexMin,
                                         TIndex indexSentinel,
                                         long stepSize,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(DynamicIndexPartitionFactory({indexMin, indexSentinel, stepSize}, _pool.size()),
                  rFunction,
                  true);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexSentinel,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(DynamicIndexPartitionFactory({0, indexSentinel, 1}, _pool.size()),
                  rFunction,
                  true);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(Ref<TContainer> rContainer,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(DynamicIndexPartitionFactory({0, rContainer.size(), 1}, _pool.size()),
                  rContainer,
                  rFunction);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(Ref<const TContainer> rContainer,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(DynamicIndexPartitionFactory({0, rContainer.size(), 1}, _pool.size()),
                  rContainer,
                  rFunction);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Iterator TIterator, class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIterator itBegin,
                                         TIterator itEnd,
                                         Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(DynamicIndexPartitionFactory({0, static_cast<std::size_t>(std::distance(itBegin, itEnd)), 1}, _pool.size()),
                  itBegin,
                  rFunction);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
Ref<const typename ParallelFor<TIndex,TStorage>::Pool>
ParallelFor<TIndex,TStorage>::getPool() const noexcept
{
    return _pool;
}


template <concepts::Integer TIndex, class TStorage>
Ref<typename ParallelFor<TIndex,TStorage>::Pool>
ParallelFor<TIndex,TStorage>::getPool() noexcept
{
    return _pool;
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                                      Ref<const TFunction> rFunction,
                                      bool barrier)
{
    CIE_BEGIN_EXCEPTION_TRACING

    const Size partitionCount = rIndexPartitionFactory.size();
    if (!partitionCount) return;

    // Schedule jobs
    if (barrier) {
        // The function is captured by reference because this function imposes
        // a barrier before it exits, so the function is guaranteed to stay alive
        // for the duration of the loop.
        for (TIndex iPartition=0; iPartition<partitionCount; ++iPartition) {
            const auto partition = rIndexPartitionFactory.getPartition(iPartition);
            _pool.queueTLSJob(
                [partition, &rFunction](auto&... rArgs) -> void {
                    for (TIndex i=partition.begin; i<partition.end; i+=partition.step) {
                        rFunction(i, rArgs...);
                    }
                }
            ); // queueJob
        } // for iPartition in range(partitionCount)
        _pool.barrier();
    } else {
        // The function must be copied into the lambda because the jobs are
        // not guaranteed to finish before exiting this function.
        for (TIndex iPartition=0; iPartition<partitionCount; ++iPartition) {
            const auto partition = rIndexPartitionFactory.getPartition(iPartition);
            _pool.queueTLSJob(
                [partition, rFunction](auto&... rArgs) -> void {
                    for (TIndex i=partition.begin; i<partition.end; i+=partition.step) {
                        rFunction(i, rArgs...);
                    }
                }
            ); // queueJob
        } // for iPartition in range(partitionCount)
    }

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                                      Ref<TContainer> rContainer,
                                      Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto partitionCount = rIndexPartitionFactory.size();
    if (!partitionCount) return;

    // Schedule jobs
    // The function and container are captured by reference because this function imposes
    // a barrier before it exits, so the captured objects are guaranteed to stay alive
    // for the duration of the loop.
    for (TIndex iPartition=0; iPartition<partitionCount; ++iPartition) {
        const auto partition = rIndexPartitionFactory.getPartition(iPartition);
        _pool.queueTLSJob(
            [partition, &rFunction, &rContainer](auto&... rArgs) -> void {
                for (TIndex i=partition.begin; i<partition.end; i+=partition.step) {
                    rFunction(rContainer[i], rArgs...);
                }
            }
        ); // queueJob
    } // for iPartition in range(partitionCount)

    _pool.barrier();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                                      Ref<const TContainer> rContainer,
                                      Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto partitionCount = rIndexPartitionFactory.size();
    if (!partitionCount) return;

    // Schedule jobs
    // The function and container are captured by reference because this function imposes
    // a barrier before it exits, so the captured objects are guaranteed to stay alive
    // for the duration of the loop.
    for (TIndex iPartition=0; iPartition<partitionCount; ++iPartition) {
        const auto partition = rIndexPartitionFactory.getPartition(iPartition);
        _pool.queueTLSJob(
            [partition, &rFunction, &rContainer](auto&... rArgs) -> void {
                for (TIndex i=partition.begin; i<partition.end; i+=partition.step) {
                    rFunction(rContainer[i], rArgs...);
                }
            }
        ); // queueJob
    } // for iPartition in range(partitionCount)

    _pool.barrier();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Iterator TIterator, class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                                      TIterator itBegin,
                                      Ref<const TFunction> rFunction)
{
    CIE_BEGIN_EXCEPTION_TRACING
    const auto partitionCount = rIndexPartitionFactory.size();
    if (!partitionCount) return;

    // Schedule jobs
    for (TIndex iPartition=0; iPartition<partitionCount; ++iPartition) {
        const auto partition = rIndexPartitionFactory.getPartition(iPartition);
        _pool.queueTLSJob(
            [partition, itBegin, &rFunction] (auto&... rArgs) -> void {
                for (TIndex i=partition.begin; i<partition.end; i+= partition.step) {
                    auto it = itBegin;
                    std::ranges::advance(it, i);
                    rFunction(*it, rArgs...);
                }
            }
        );
    }

    _pool.barrier();
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::mp


#endif
