#ifndef CIE_UTILS_PARALLEL_FOR_IMPL_HPP
#define CIE_UTILS_PARALLEL_FOR_IMPL_HPP

// --- Utility Includes ---
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
ParallelFor<TIndex,TStorage>::ParallelFor(Ref<ThreadPoolBase> r_pool, TArgs&&... r_args) :
    _pool(r_pool, TStorage(std::forward<TArgs>(r_args)...))
{
}


template <concepts::Integer TIndex, class TStorage>
template <class ...TArgs>
ParallelFor<TIndex,ThreadStorage<typename std::remove_reference_t<TArgs>...>>
ParallelFor<TIndex,TStorage>::firstPrivate(TArgs&&... r_args)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return ParallelFor<
        TIndex,
        ThreadStorage<typename std::remove_reference_t<TArgs>...>
    >(_pool, std::forward<TArgs>(r_args)...);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexMin,
                                         TIndex indexMax,
                                         TIndex stepSize,
                                         Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(this->makeIndexPartition(indexMin, indexMax, stepSize),
                  stepSize,
                  r_function);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexMax,
                                         Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(this->makeIndexPartition(0, indexMax, 1),
                  1,
                  r_function);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(Ref<TContainer> r_container,
                                         Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(this->makeIndexPartition(0,r_container.size(),1),
                  1,
                  r_container,
                  r_function);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(Ref<const TContainer> r_container,
                                         Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->execute(this->makeIndexPartition(0,r_container.size(),1),
                  1,
                  r_container,
                  r_function);
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
typename ParallelFor<TIndex,TStorage>::IndexPartition
ParallelFor<TIndex,TStorage>::makeIndexPartition(TIndex indexMin,
                                                 TIndex indexMax,
                                                 TIndex stepSize)
{
    CIE_BEGIN_EXCEPTION_TRACING

    IndexPartition indexPartition;

    // Check index range
    CIE_CHECK(stepSize != 0, "0 step size in ParallelFor")
    CIE_CHECK((indexMin <= indexMax) == (0 < stepSize), "Invalid step direction in ParallelFor");

    // Handle reverse loops
    if (stepSize < 0) {
        std::swap(indexMin, indexMax);
        stepSize = -stepSize;
    }

    // Initialize index blocks
    TIndex blockSize = (indexMax - indexMin) / stepSize / _pool.size();
    blockSize += blockSize % stepSize;
    blockSize = std::max(blockSize, 1ul);

    const Size numberOfThreads = _pool.size();
    indexPartition.reserve(numberOfThreads + 1);
    TIndex tmp = indexMin;
    for (TIndex i=0; i<numberOfThreads && tmp<indexMax; ++i, tmp+=blockSize) {
        indexPartition.push_back(tmp);
    }
    indexPartition.push_back(indexMax);
    return indexPartition;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartition> r_indexPartition,
                                      TIndex stepSize,
                                      Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (r_indexPartition.empty()) return;
    CIE_OUT_OF_RANGE_CHECK(r_indexPartition.size() <= _pool.size() + 1)

    // Schedule jobs
    for (TIndex partitionIndex=0; partitionIndex<r_indexPartition.size()-1; ++partitionIndex) {
        _pool.queueTLSJob(
            [partitionIndex, stepSize, &r_indexPartition, &r_function](auto&... r_args) -> void {
                for (TIndex i=r_indexPartition[partitionIndex]; i<r_indexPartition[partitionIndex+1]; i+=stepSize) {
                    r_function(i, r_args...);
                }
            }
        ); // queueJob
    }

    _pool.barrier();
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartition> r_indexPartition,
                                      TIndex stepSize,
                                      Ref<TContainer> r_container,
                                      Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (r_indexPartition.empty()) return;
    CIE_OUT_OF_RANGE_CHECK(r_indexPartition.size() <= _pool.size() + 1)

    // Schedule jobs
    for (TIndex partitionIndex=0; partitionIndex<r_indexPartition.size()-1; ++partitionIndex) {
        _pool.queueTLSJob(
            [partitionIndex, stepSize, &r_indexPartition, &r_function, &r_container](auto&... r_args) -> void {
                for (TIndex i=r_indexPartition[partitionIndex]; i<r_indexPartition[partitionIndex+1]; i+=stepSize)
                    r_function(r_container[i], r_args...);
            }
        ); // queueJob
    }

    _pool.barrier();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer, class TFunction>
void
ParallelFor<TIndex,TStorage>::execute(Ref<const IndexPartition> r_indexPartition,
                                      TIndex stepSize,
                                      Ref<const TContainer> r_container,
                                      Ref<const TFunction> r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (r_indexPartition.empty()) return;
    CIE_OUT_OF_RANGE_CHECK(r_indexPartition.size() <= _pool.size() + 1)

    // Schedule jobs
    for (TIndex partitionIndex=0; partitionIndex<r_indexPartition.size()-1; ++partitionIndex) {
        _pool.queueTLSJob(
            [partitionIndex, stepSize, &r_indexPartition, &r_function, &r_container](auto&... r_args) -> void {
                for (TIndex i=r_indexPartition[partitionIndex]; i<r_indexPartition[partitionIndex+1]; i+=stepSize)
                    r_function(r_container[i], r_args...);
            }
        ); // queueJob
    }

    _pool.barrier();

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::mp


#endif
