#ifndef CIE_UTILS_CONCURRENCY_THREAD_POOL_BASE_IMPL_HPP
#define CIE_UTILS_CONCURRENCY_THREAD_POOL_BASE_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"


namespace cie::mp {


inline Size ThreadPoolBase::size() const noexcept
{
    return _threads.size();
}


inline Size ThreadPoolBase::numberOfJobs() const noexcept
{
    ScopedLock lock(_mutex);
    return _jobs.size();
}


inline void ThreadPoolBase::queueJob(Job&& r_job)
{
    {
        ScopedLock lock(_mutex);
        _jobs.emplace(std::move(r_job));
    }
    _idleCondition.notify_one();
}


inline void ThreadPoolBase::queueJob(const Job& r_job)
{
    {
        ScopedLock lock(_mutex);
        _jobs.push(r_job);
    }
    _idleCondition.notify_one();
}


} // namespace cie::mp


#endif
