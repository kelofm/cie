#ifndef CIE_UTILS_PARALLEL_FOR_IMPL_HPP
#define CIE_UTILS_PARALLEL_FOR_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolSingleton.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::mp {


template <concepts::Integer TIndex, class TStorage>
template <class ...TArgs>
requires (!std::is_same_v<std::tuple<typename std::decay<TArgs>::type...>,std::tuple<ParallelFor<TIndex,TStorage>>>)
ParallelFor<TIndex,TStorage>::ParallelFor(TArgs&&... r_args) :
    _pool(ThreadPoolSingleton::firstPrivate(std::forward<TArgs>(r_args)...))
{
}


template <concepts::Integer TIndex, class TStorage>
template <class ...TArgs>
inline ParallelFor<TIndex,ThreadStorage<typename std::remove_reference<TArgs>::type...>>
ParallelFor<TIndex,TStorage>::firstPrivate(TArgs&&... r_args)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return ParallelFor<
        TIndex,
        ThreadStorage<typename std::remove_reference<TArgs>::type...>
    >(std::forward<TArgs>(r_args)...);

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
inline void
ParallelFor<TIndex,TStorage>::setPool(RightRef<Pool> r_pool) noexcept
{
    _pool = std::move(r_pool);
}


template <concepts::Integer TIndex, class TStorage>
inline void
ParallelFor<TIndex,TStorage>::setPool(Ref<Pool> r_pool)
{
    _pool = r_pool;
}


template <concepts::Integer TIndex, class TStorage>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexMin,
                                         TIndex indexMax,
                                         TIndex stepSize,
                                         IndexLoopFunction&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(indexMin,indexMax,stepSize)),
        stepSize,
        std::forward<IndexLoopFunction>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TIndex indexMax,
                                         IndexLoopFunction&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(0, indexMax, 1)),
        1,
        std::forward<IndexLoopFunction>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TContainer& r_container,
                                         typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(0,r_container.size(),1)),
        1,
        r_container,
        std::forward<typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(const TContainer& r_container,
                                         typename TStorage::template ConstObjectLoopFunction<typename TContainer::value_type>&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(0,r_container.size(),1)),
        1,
        r_container,
        std::forward<typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(TContainer& r_container,
                                         typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(0,r_container.size(),1)),
        1,
        r_container,
        std::forward<typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer>
inline ParallelFor<TIndex,TStorage>&
ParallelFor<TIndex,TStorage>::operator()(const TContainer& r_container,
                                         typename TStorage::template ConstObjectLoopFunction<typename TContainer::value_type>& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->execute(
        std::move(this->makeIndexPartition(0,r_container.size(),1)),
        1,
        r_container,
        std::forward<typename TStorage::template ObjectLoopFunction<typename TContainer::value_type>>(r_function)
    );

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
inline Ref<const typename ParallelFor<TIndex,TStorage>::Pool>
ParallelFor<TIndex,TStorage>::getPool() const noexcept
{
    return _pool;
}


template <concepts::Integer TIndex, class TStorage>
inline Ref<typename ParallelFor<TIndex,TStorage>::Pool>
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

    // Check index range
    CIE_OUT_OF_RANGE_CHECK(stepSize != 0)
    CIE_OUT_OF_RANGE_CHECK((indexMin <= indexMax) == (0 < stepSize));

    // Handle reverse loops
    if (stepSize < 0) {
        std::swap(indexMin, indexMax);
        stepSize = -stepSize;
    }

    typename ParallelFor<TIndex,TStorage>::IndexPartition indexPartition;

    // Initialize index blocks
    TIndex blockSize = stepSize * (indexMax - indexMin) / (stepSize * _pool.size() + 1);
    blockSize = std::max(blockSize, 1UL);

    TIndex tmp = indexMin;

    const Size numberOfThreads = _pool.size();
    indexPartition.reserve(numberOfThreads + 1);
    for (TIndex i=0; i<numberOfThreads && tmp<indexMax; ++i, tmp+=blockSize)
        indexPartition.push_back(tmp);

    indexPartition.push_back(indexMax);

    return indexPartition;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
void
ParallelFor<TIndex,TStorage>::execute(const typename ParallelFor<TIndex,TStorage>::IndexPartition& r_indexPartition,
                                      TIndex stepSize,
                                      typename ParallelFor<TIndex,TStorage>::IndexLoopFunction&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_indexPartition.size() <= _pool.size() + 1)

    // Schedule jobs
    for (TIndex partitionIndex=0; partitionIndex<r_indexPartition.size()-1; ++partitionIndex)
    {
        _pool.queueTLSJob(
            [partitionIndex, stepSize, r_indexPartition, r_function](auto&... r_args) -> void {
                for (TIndex i=r_indexPartition[partitionIndex]; i<r_indexPartition[partitionIndex+1]; i+=stepSize)
                    r_function(i, r_args...);
            }
        ); // queueJob
    }

    _pool.barrier();

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, class TStorage>
template <concepts::Container TContainer>
void
ParallelFor<TIndex,TStorage>::execute(const typename ParallelFor<TIndex,TStorage>::IndexPartition& r_indexPartition,
                                      TIndex stepSize,
                                      TContainer& r_container,
                                      typename ParallelFor<TIndex,TStorage>::ObjectLoopFunction<typename TContainer::value_type>&& r_function)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_indexPartition.size() <= _pool.size() + 1)

    // Schedule jobs
    for (TIndex partitionIndex=0; partitionIndex<r_indexPartition.size()-1; ++partitionIndex) {
        _pool.queueTLSJob(
            [partitionIndex, stepSize, r_indexPartition, r_function, &r_container](auto&... r_args) -> void {
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
