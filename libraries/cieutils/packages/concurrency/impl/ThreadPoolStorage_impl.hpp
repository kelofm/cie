#ifndef CIE_UTILS_THREAD_POOL_STORAGE_IMPL_HPP
#define CIE_UTILS_THREAD_POOL_STORAGE_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolStorage.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::mp {


template <class TStorage>
ThreadPoolStorage<TStorage>::ThreadPoolStorage(Ref<ThreadPoolBase> r_pool)
    : ThreadPoolStorage(r_pool, TStorage())
{
}


template <class TStorage>
ThreadPoolStorage<TStorage>::ThreadPoolStorage(Ref<ThreadPoolBase> r_pool, Ref<const TStorage> r_storage)
    : _r_pool(r_pool),
      _storages()
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto thisThreadID = std::this_thread::get_id();

    for (const auto& r_thread : _r_pool.get().threads())
        _storages.emplace(r_thread.get_id(), r_storage);

    // This object must be constructed from the master thread!
    CIE_CHECK(_storages.find(thisThreadID) == _storages.end(),
              "ThreadPoolStorage constructed from non-master thread " << thisThreadID)

    _storages.emplace(thisThreadID, r_storage);

    CIE_END_EXCEPTION_TRACING
}


template <class TStorage>
inline Ref<TStorage> ThreadPoolStorage<TStorage>::getLocalStorage()
{
    return _storages[std::this_thread::get_id()];
}


template <class TStorage>
inline Ref<const TStorage> ThreadPoolStorage<TStorage>::getLocalStorage() const
{
    return _storages[std::this_thread::get_id()];
}


template <class TStorage>
void ThreadPoolStorage<TStorage>::setLocalStorage(RightRef<TStorage> r_storage)
{
    _storages[std::this_thread::get_id()] = std::move(r_storage);
}


template <class TStorage>
void ThreadPoolStorage<TStorage>::setLocalStorage(const TStorage& r_storage)
{
    _storages[std::this_thread::get_id()] = r_storage;
}


} // namespace cie::mp


#endif