#ifndef CIE_UTILS_CONCURRENCY_THREAD_POOL_BASE_HPP
#define CIE_UTILS_CONCURRENCY_THREAD_POOL_BASE_HPP

// --- Internal Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- STL Includes ---
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <queue>
#include <functional>
#include <memory>
#include <atomic>
#include <iosfwd>


namespace cie::mp {


class ThreadPoolBase final
{
public:
    using Job = std::function<void()>;

    CIE_DEFINE_CLASS_POINTERS(ThreadPoolBase)

private:
    using Thread = std::thread;

    using ThreadContainer = std::vector<Thread>;

    using JobContainer = std::queue<Job>;

    using Mutex = std::mutex;

    using Lock = std::unique_lock<Mutex>;

    using ScopedLock = std::scoped_lock<Mutex>;

public:
    /// @brief Create a thread pool with 'size' threads
    /**
     *  @note @a size must be positive and not greater than the maximum number of threads
     *         on the system. The number of threads is capped at the supported max.
     */
    ThreadPoolBase(Size size);

    /// Create a thread pool with the maximum number of threads the system supports
    ThreadPoolBase();

    ~ThreadPoolBase();

    /// @brief Get the maximum number of threads supported on the system.
    static Size maxNumberOfThreads() noexcept;

    /// @brief Queue a job for execution
    void queueJob(Job&& r_job);

    /// @brief Queue a job for execution
    void queueJob(const Job& r_job);

    /// @brief Number of threads in the pool
    Size size() const noexcept;

    /// @brief Number of jobs in the queue.
    Size numberOfJobs() const noexcept;

    /// @brief Block execution until all threads reach this barrier.
    void barrier();

    /// @brief Get the index of the caller's thread in the internal container.
    Size threadIndex() const;

    /// @brief Get immutable access to the internal container of threads.
    Ref<const ThreadContainer> threads() const;

private:
    ThreadPoolBase(ThreadPoolBase&& r_rhs) = delete;

    ThreadPoolBase(const ThreadPoolBase& r_rhs) = delete;

    ThreadPoolBase& operator=(ThreadPoolBase&& r_rhs) = delete;

    ThreadPoolBase& operator=(const ThreadPoolBase& r_rhs) = delete;

    /** @brief Looping function running on all threads
     *  @details Upon notification, sets the lock and checks the job queue,
     *           strips and executes one if there are any.
     */
    void jobScheduler();

    /// @brief Clear all scheduled jobs and join all threads after they finished executing their current jobs.
    void terminate();

private:
    std::atomic<bool> _terminate;

    ThreadContainer _threads;

    JobContainer _jobs;

    mutable Mutex _mutex;

    std::condition_variable _idleCondition;

    std::condition_variable _barrierCondition;

    std::condition_variable _masterCondition;

    std::atomic<Size> _numberOfActiveThreads;

    std::atomic<bool> _barrier;
};


Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const ThreadPoolBase> r_pool);


} // namespace cie::mp

#include "packages/concurrency/impl/ThreadPoolBase_impl.hpp"

#endif
