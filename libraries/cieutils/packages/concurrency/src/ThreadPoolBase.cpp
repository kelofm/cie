// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <algorithm>
#include <cstdlib>
#include <stdlib.h>


namespace cie::mp {


ThreadPoolBase::ThreadPoolBase(Size size) :
    _terminate(false)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (size == 0)
        CIE_THROW(Exception, "Cannot create thread pool of size 0!")

    // Get number of threads available in the machine
    Size maxNumberOfThreads = ThreadPoolBase::maxNumberOfThreads();

    // Limit max number of threads if the environment defines OMP_NUM_THREADS
    const auto p_ompNumThreads = std::getenv("OMP_NUM_THREADS");
    if (p_ompNumThreads) {
        const Size ompNumThreads = std::strtoul(p_ompNumThreads, NULL, 10);
        if (ompNumThreads < 1)
            CIE_THROW(Exception, "Invalid OMP_NUM_THREADS (" + std::to_string(ompNumThreads) + ")")

        if (ompNumThreads < maxNumberOfThreads)
            maxNumberOfThreads = ompNumThreads;
    }

    if (maxNumberOfThreads == 0)
        maxNumberOfThreads = 1;

    if (maxNumberOfThreads < size)
        size = maxNumberOfThreads;

    _barrier = false;

    // Initialize threads
    _threads.reserve(size);

    for (_numberOfActiveThreads=0; _numberOfActiveThreads<size; ++_numberOfActiveThreads)
        _threads.emplace_back(&ThreadPoolBase::jobScheduler, this);

    CIE_END_EXCEPTION_TRACING
}


ThreadPoolBase::ThreadPoolBase() :
    ThreadPoolBase(ThreadPoolBase::maxNumberOfThreads())
{
}


ThreadPoolBase::~ThreadPoolBase()
{
    terminate();
}


Size ThreadPoolBase::maxNumberOfThreads() noexcept
{
    return std::thread::hardware_concurrency();
}


void ThreadPoolBase::barrier()
{
    CIE_BEGIN_EXCEPTION_TRACING

    Lock lock(_mutex);

    CIE_CHECK(
        !_barrier,
        "barrier was called on an ongoing barrier!"
    )

    _masterCondition.wait(
        lock,
        [this]{return _numberOfActiveThreads == _threads.size();}
    );

    _barrier = true;

    _idleCondition.notify_all();

    _masterCondition.wait(
        lock,
        [this]{return (_terminate) || (_numberOfActiveThreads == 0);}
    );

    _barrier = false;
    _barrierCondition.notify_all();

    CIE_END_EXCEPTION_TRACING
}


void ThreadPoolBase::terminate()
{
    {
        ScopedLock lock(_mutex);
        _terminate = true;
        _idleCondition.notify_all();
        _barrierCondition.notify_all();
    }

    for (auto& r_thread : _threads)
        r_thread.join();

    _threads.clear();
}


void ThreadPoolBase::jobScheduler()
{
    while (true)
    {
        ThreadPoolBase::Job job = nullptr;

        {
            Lock lock(_mutex);

            _idleCondition.wait(
                lock,
                [this]{return !_jobs.empty() || _terminate || _barrier;}
            );

            // Got a job? -> strip it!
            if (!_jobs.empty())
            {
                job = _jobs.front();
                _jobs.pop();
            }
        }

        // Have a stripped job? -> execute it!
        if (job)
            job();

        // Don't have a stripped job? -> terminate or set blocked for barrier
        else
        {
            Lock lock(_mutex);

            // Terminate thread
            if (_terminate)
                break;

            // Enqueue at barrier
            else if (_barrier)
            {
                --_numberOfActiveThreads;

                _masterCondition.notify_all();

                _barrierCondition.wait(
                    lock,
                    [this]{ return (!_barrier) || _terminate; }
                );

                ++_numberOfActiveThreads;
                _masterCondition.notify_all();
            }
            else
                CIE_THROW(Exception, "Thread has no job, no barrier, and is not terminated")
        } // job queue empty
    }
}


Size ThreadPoolBase::threadIndex() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    const auto id = std::this_thread::get_id();
    const auto it_end = _threads.cend();
    const auto it_thread = std::find_if(_threads.begin(),
                                        it_end,
                                        [id](const Thread& r_thread)
                                            {return id == r_thread.get_id();});

    CIE_CHECK(it_thread != it_end,
              "Thread " << id << " is not part of the thread pool")

    return std::distance(_threads.begin(), it_thread);

    CIE_END_EXCEPTION_TRACING
}


Ref<const ThreadPoolBase::ThreadContainer> ThreadPoolBase::threads() const
{
    return _threads;
}


} // namespace cie::mp