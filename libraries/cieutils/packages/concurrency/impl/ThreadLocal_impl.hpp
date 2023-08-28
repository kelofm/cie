#ifndef CIE_UTILS_THREAD_LOCAL_IMPL_HPP
#define CIE_UTILS_THREAD_LOCAL_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadLocal.hpp"
#include "packages/concurrency/inc/ThreadUtils.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include <thread>


namespace cie::mp {


template <class ...TStored>
ThreadLocal<TStored...>::ThreadLocal()
{
    for (std::thread::id id : getThreadIDs()) {
        _map.try_emplace(id);
    } // for thread.id in threads

    if constexpr (0 < std::tuple_size_v<Tuple>) {
        auto it_locals = _map.find(std::this_thread::get_id());
        CIE_OUT_OF_RANGE_CHECK(it_locals != _map.end())
        this->populateLocals<0>(it_locals);
    }
}


template <class ...TStored>
ThreadLocal<TStored...>::ThreadLocal(const TStored&... r_items)
    : ThreadLocal()
{
    CIE_BEGIN_EXCEPTION_TRACING
    const auto it_end = _map.end();
    for (auto it_pair=_map.begin(); it_pair!=it_end; it_pair++) {
        it_pair.value() = Tuple(r_items ...);
    } // for pair in map

    if constexpr (0 < std::tuple_size_v<Tuple>) {
        auto it_locals = _map.find(std::this_thread::get_id());
        CIE_OUT_OF_RANGE_CHECK(it_locals != _map.end())
        this->populateLocals<0>(it_locals);
    }
    CIE_END_EXCEPTION_TRACING
}


template <class ...TStored>
template <std::size_t Index>
Ref<std::tuple_element_t<Index,typename ThreadLocal<TStored...>::Tuple>>
inline ThreadLocal<TStored...>::get()
{
    return *std::get<Index>(_locals);
}


template <class ...TStored>
template <std::size_t Index>
Ref<const std::tuple_element_t<Index,typename ThreadLocal<TStored...>::Tuple>>
inline ThreadLocal<TStored...>::get() const
{
    return *std::get<Index>(_locals);
}


template <class ...TStored>
template <unsigned Index>
void ThreadLocal<TStored...>::populateLocals(typename Map::iterator it_locals)
{
    std::get<Index>(_locals) = &std::get<Index>(it_locals.value());
    if constexpr (Index + 1 < std::tuple_size_v<Tuple>) {
        this->populateLocals<Index+1>(it_locals);
    }
}


} // namespace cie::mp


#endif
