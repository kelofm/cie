#ifndef CIE_UTILS_THREAD_LOCAL_HPP
#define CIE_UTILS_THREAD_LOCAL_HPP

// --- External Includes ---
#include "tsl/robin_map.h"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <tuple>
#include <thread>


namespace cie::mp {


/// @addtogroup cieutils
template <class ...TStored>
class ThreadLocal
{
public:
    using Tuple = std::tuple<TStored...>;

    using Map = tsl::robin_map<std::thread::id,Tuple>;

public:
    ThreadLocal();

    ThreadLocal(const TStored&... r_items);

    template <std::size_t Index>
    Ref<std::tuple_element_t<Index,Tuple>> get();

    template <std::size_t Index>
    Ref<const std::tuple_element_t<Index,Tuple>> get() const;

private:
    /// @brief Fill the local pointer storage.
    template <unsigned Index>
    void populateLocals(typename Map::iterator it_locals);

private:
    /// @details Store pointers to the local objects to avoid having to
    ///          go through the hash map each time an object is requested.
    std::tuple<Ptr<TStored>...> _locals;

    Map _map;
}; // class ThreadLocal


} // namespace cie::mp


#endif

#include "packages/concurrency/impl/ThreadLocal_impl.hpp"
