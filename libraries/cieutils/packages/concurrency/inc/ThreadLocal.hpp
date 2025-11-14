#pragma once

// --- External Includes ---
#include "tsl/robin_map.h"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <tuple>
#include <thread>


namespace cie::mp {


/// @ingroup cieutils
template <class ...TStored>
class ThreadLocal
{
public:
    using Tuple = std::tuple<TStored...>;

    using Map = tsl::robin_map<std::thread::id,Tuple>;

public:
    ThreadLocal(Ref<const ThreadPoolBase> rPool);

    ThreadLocal(Ref<const ThreadPoolBase> rPool,
                const TStored&... rItems);

    template <std::size_t Index>
    Ref<std::tuple_element_t<Index,Tuple>> get();

    template <std::size_t Index>
    Ref<const std::tuple_element_t<Index,Tuple>> get() const;

    Ref<const ThreadPoolBase> threadPool() const noexcept;

private:
    Ptr<const ThreadPoolBase> _pThreadPool;

    Map _map;
}; // class ThreadLocal


} // namespace cie::mp

#include "packages/concurrency/impl/ThreadLocal_impl.hpp"
