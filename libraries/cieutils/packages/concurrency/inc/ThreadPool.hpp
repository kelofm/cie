#ifndef CIE_UTILS_THREAD_POOL_HPP
#define CIE_UTILS_THREAD_POOL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/concurrency/inc/ThreadPoolStorage.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <utility> // std::reference_wrapper
#include <iosfwd>


namespace cie::mp {


/// @addtogroup cieutils
template <class TStorage = ThreadStorage<>>
class ThreadPool final
{
public:
    /// @brief Standard job type queued in @ref ThreadPoolBase.
    using Job = ThreadPoolBase::Job;

    /// @brief Job type accepting arguments from the thread local storage @a TStorage.
    using TLSJob = typename TStorage::template Function<void>;

    /// @brief Thread local storage type.
    using Storage = TStorage;

    CIE_DEFINE_CLASS_POINTERS(ThreadPool)

public:
    ThreadPool(Ref<ThreadPoolBase> r_pool);

    ThreadPool(const TStorage& r_storage);

    ThreadPool(Ref<ThreadPoolBase> r_pool, RightRef<TStorage> r_storage);

    ThreadPool(Ref<ThreadPoolBase> r_pool, Ref<const TStorage> r_storage);

    ThreadPool(ThreadPool&& r_rhs) noexcept = default;

    ThreadPool(const ThreadPool& r_rhs) = default;

    ThreadPool& operator=(ThreadPool&& r_rhs) noexcept = default;

    ThreadPool& operator=(const ThreadPool& r_rhs);

    template <class ...Ts>
    ThreadPool<ThreadStorage<Ts...>>
    firstPrivate(Ts&&... r_storage);

    /// @brief Get the maximum number of theads supported on the system.
    static Size maxNumberOfThreads() noexcept;

    /// @brief Queue a job using variables stored in thread local storage.
    void queueTLSJob(RightRef<TLSJob> r_job);

    /// @brief Queue a job using variables stored in thread local storage.
    void queueTLSJob(Ref<const TLSJob> r_job);

    /// @brief Directly queue a job in the wrapped @ref ThreadPoolBase.
    void queueJob(RightRef<Job> r_job);

    /// @brief Directly queue a job in the wrapped @ref ThreadPoolBase.
    void queueJob(Ref<const Job> r_job);

    /// @brief Number of threads in the wrapped @ref ThreadPoolBase.
    Size size() const noexcept;

    /// @brief Number of jobs queued in the wrapped @ref ThreadPoolBase.
    Size numberOfJobs() const noexcept;

    /// @brief Issue a barrier on the wrapped @ref ThreadPoolBase.
    void barrier();

    /// @brief Provides mutable access to the thread local storage of all threads.
    Ref<ThreadPoolStorage<TStorage>> getStorage() noexcept;

    /// @brief Provides immutable access to the thread local storage of all threads.
    Ref<const ThreadPoolStorage<TStorage>> getStorage() const noexcept;

    operator const ThreadPoolBase& () const noexcept;

    operator ThreadPoolBase& () noexcept;

    template <class TS>
    friend Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const ThreadPool<TS>> r_pool);

private:
    std::reference_wrapper<ThreadPoolBase> _r_pool;

    ThreadPoolStorage<TStorage> _storage;
}; // class ThreadPool


} // namespace cie::mp


namespace cie::concepts {

/// @addtogroup cieutils
template <class T, class ...TStored>
concept ThreadPool
= requires (typename std::decay<T>::type instance)
{
    typename T::Job;

    typename T::TLSJob;

    typename T::Storage;

    {instance.queueTLSJob(typename T::TLSJob())};

    {instance.queueJob(typename T::Job())};

    {instance.barrier()};
} ;//&& std::is_same_v<typename T::Storage,mp::ThreadStorage<TStored...>>; // concept ThreadPool


} // namespace cie::concepts

#include "packages/concurrency/impl/ThreadPool_impl.hpp"

#endif
