// --- Internal Includes ---
#include "packages/concurrency/inc/ThreadPoolSingleton.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <thread>
#include <algorithm>


namespace cie::mp {


std::optional<ThreadPoolBase::UniquePointer> ThreadPoolSingleton::_p_pool;


Ref<ThreadPoolBase> ThreadPoolSingleton::getBase()
{
    if (!_p_pool)
        _p_pool.emplace(new ThreadPoolBase);

    return *_p_pool.value();
}


void ThreadPoolSingleton::shutDown()
{
    if (_p_pool) [[likely]] {
         _p_pool.reset();
    }
}


ThreadPool<> ThreadPoolSingleton::get()
{
    return ThreadPool<>(ThreadPoolSingleton::getBase());
}


} // namespace cie::mp
