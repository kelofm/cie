#ifndef CIE_UTILS_THREAD_POOL_IMPL_HPP
#define CIE_UTILS_THREAD_POOL_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPool.hpp"


namespace cie::mp {


class ThreadPoolSingleton;


template <class TStorage>
ThreadPool<TStorage>::ThreadPool(Ref<ThreadPoolBase> r_pool)
    : ThreadPool(r_pool, TStorage())
{
}


template <class TStorage>
ThreadPool<TStorage>::ThreadPool(Ref<ThreadPoolBase> r_pool, RightRef<TStorage> r_storage)
    : _r_pool(r_pool),
      _storage(r_pool, std::move(r_storage))
{
}


template <class TStorage>
ThreadPool<TStorage>::ThreadPool(Ref<ThreadPoolBase> r_pool, Ref<const TStorage> r_storage)
    : _r_pool(r_pool),
      _storage(r_pool, r_storage)
{
}


template <class TStorage>
ThreadPool<TStorage>&
ThreadPool<TStorage>::operator=(const ThreadPool& r_rhs)
{
    _r_pool = r_rhs._r_pool.get();
    _storage = r_rhs._storage;
    return *this;
}



template <class TStorage>
inline Size
ThreadPool<TStorage>::maxNumberOfThreads() noexcept
{
    return ThreadPoolBase::maxNumberOfThreads();
}


template <class TStorage>
inline void
ThreadPool<TStorage>::queueTLSJob(RightRef<TLSJob> r_job)
{
    _r_pool.get().queueJob(
        [this, r_job = std::move(r_job)]()
        {
            std::apply(
                [r_job](auto&... r_args) -> void {r_job(r_args...);},
                _storage.getLocalStorage().getRefs()
            );
        }
    );
}


template <class TStorage>
inline void
ThreadPool<TStorage>::queueTLSJob(Ref<const TLSJob> r_job)
{
    _r_pool.get().queueJob(
        [this, r_job]()
        {
            std::apply(
                [&r_job](auto&... r_args) -> void {r_job(r_args...);},
                _storage.getLocalStorage().getRefs()
            );
        }
    );
}


template <class TStorage>
inline void ThreadPool<TStorage>::queueJob(RightRef<Job> r_job)
{
    _r_pool.get().queueJob(std::move(r_job));
}


template <class TStorage>
inline void ThreadPool<TStorage>::queueJob(Ref<const Job> r_job)
{
    _r_pool.get().queueJob(std::move(r_job));
}


template <class TStorage>
inline Size
ThreadPool<TStorage>::size() const noexcept
{
    return _r_pool.get().size();
}


template <class TStorage>
inline Size
ThreadPool<TStorage>::numberOfJobs() const noexcept
{
    return _r_pool.get().numberOfJobs();
}


template <class TStorage>
void
ThreadPool<TStorage>::barrier()
{
    _r_pool.get().barrier();
}


template <class TStorage>
inline Ref<ThreadPoolStorage<TStorage>>
ThreadPool<TStorage>::getStorage() noexcept
{
    return _storage;
}


template <class TStorage>
inline Ref<const ThreadPoolStorage<TStorage>>
ThreadPool<TStorage>::getStorage() const noexcept
{
    return _storage;
}



template <class TStorage>
Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const ThreadPool<TStorage>> r_pool)
{
    return r_stream << r_pool._r_pool;
}


} // namespace cie::fem


#endif