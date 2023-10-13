#ifndef CIE_UTILS_THREAD_POOL_SINGLETON_HPP
#define CIE_UTILS_THREAD_POOL_SINGLETON_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"

// --- STL Includes ---
#include <optional>


namespace cie::mp {


class ThreadPoolSingleton
{
public:
    [[nodiscard]] static Ref<ThreadPoolBase> getBase();

    static void clear();

    [[nodiscard]] static ThreadPool<> get();

    template <class ...TArgs>
    [[nodiscard]] static ThreadPool<ThreadStorage<typename std::remove_reference<TArgs>::type...>> firstPrivate(TArgs&&... r_args);

private:
    static std::optional<ThreadPoolBase::UniquePointer> _p_pool;
}; // class ThreadPoolSingleton


} // namespace cie::mp

#include "packages/concurrency/impl/ThreadPoolSingleton_impl.hpp"

#endif