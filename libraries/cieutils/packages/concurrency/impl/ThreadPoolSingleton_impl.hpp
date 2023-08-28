#ifndef CIE_UTILS_THREAD_POOL_SINGLETON_IMPL_HPP
#define CIE_UTILS_THREAD_POOL_SINGLETON_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolSingleton.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::mp {


template <class ...TArgs>
ThreadPool<ThreadStorage<typename std::remove_reference<TArgs>::type...>>
ThreadPoolSingleton::firstPrivate(TArgs&&... r_args)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using StorageType = ThreadStorage<typename std::remove_reference<TArgs>::type...>;
    auto& r_base = ThreadPoolSingleton::getBase();

    return ThreadPool<StorageType>(
            r_base,
            StorageType(std::forward<TArgs>(r_args)...)
    );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::mp


#endif
