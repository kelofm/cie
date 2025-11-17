#pragma once

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadLocal.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include <packages/exceptions/inc/exception.hpp>

// --- STL Includes ---
#include <thread>


namespace cie::mp {


template <class ...TStored>
ThreadLocal<TStored...>::ThreadLocal(Ref<const ThreadPoolBase> rPool)
    : _pThreadPool(&rPool),
      _map()
{
    CIE_BEGIN_EXCEPTION_TRACING
    _map.try_emplace(std::this_thread::get_id());
    for (Ref<const std::thread> rThread : rPool.threads()) {
        _map.try_emplace(rThread.get_id());
    } // for thread.id in threads
    CIE_END_EXCEPTION_TRACING
}


template <class ...TStored>
ThreadLocal<TStored...>::ThreadLocal(Ref<const ThreadPoolBase> rPool,
                                     const TStored&... rItems)
    : ThreadLocal(rPool)
{
    CIE_BEGIN_EXCEPTION_TRACING
    _map[std::this_thread::get_id()] = Tuple(rItems...);
    for (Ref<const std::thread> rThread : rPool.threads()) {
        _map[rThread.get_id()] = Tuple(rItems...);
    }
    CIE_END_EXCEPTION_TRACING
}


template <class ...TStored>
template <std::size_t Index>
Ref<std::tuple_element_t<Index,typename ThreadLocal<TStored...>::Tuple>>
ThreadLocal<TStored...>::get()
{
    const auto it = _map.find(std::this_thread::get_id());
    CIE_OUT_OF_RANGE_CHECK(it != _map.end())
    return std::get<Index>(it.value());
}


template <class ...TStored>
template <std::size_t Index>
Ref<const std::tuple_element_t<Index,typename ThreadLocal<TStored...>::Tuple>>
ThreadLocal<TStored...>::get() const
{
    const auto it = _map.find(std::this_thread::get_id());
    CIE_OUT_OF_RANGE_CHECK(it != _map.end())
    return std::get<Index>(it.value());
}


template <class ...TStored>
Ref<const ThreadPoolBase> ThreadLocal<TStored...>::threadPool() const noexcept
{
    return *_pThreadPool;
}


} // namespace cie::mp
