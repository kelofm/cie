#ifndef CIE_UTILS_CONCURRENCY_THREAD_STORAGE_HPP
#define CIE_UTILS_CONCURRENCY_THREAD_STORAGE_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <tuple>
#include <functional>


namespace cie::mp {


template <class ...TArgs>
class ThreadStorage
{
public:
    using StorageType = std::tuple<TArgs...>;

    using RefStorage = std::tuple<TArgs&...>;

    using ConstRefStorage = std::tuple<const TArgs&...>;

    template <class TReturn>
    using Function = std::function<TReturn(TArgs&...)>;

    template <class TReturn>
    using ConstFunction = std::function<TReturn(const TArgs&...)>;

    template <concepts::Integer TIndex>
    using IndexLoopFunction = std::function<void(TIndex,TArgs&...)>;

    template <class TIndex>
    using ConstIndexLoopFunction = std::function<void(TIndex,const TArgs&...)>;

    template <class TValue>
    using ObjectLoopFunction = std::function<void(TValue&,TArgs&...)>;

    template <class TValue>
    using ConstObjectLoopFunction = std::function<void(const TValue&, TArgs&...)>;

public:
    ThreadStorage(ThreadStorage&& r_rhs) = default;

    ThreadStorage(const ThreadStorage& r_rhs) = default;

    template <class ...Ts>
    ThreadStorage(Ts&&... r_args) :
        _values(std::forward<Ts>(r_args)...)
    {}

    ThreadStorage& operator=(ThreadStorage&& r_rhs) = default;

    ThreadStorage& operator=(const ThreadStorage& r_rhs) = default;

    RefStorage getRefs()
    {return std::make_from_tuple<RefStorage>(_values);}

    ConstRefStorage getRefs() const
    {return std::make_from_tuple<ConstRefStorage>(_values);}

    template <Size Index>
    typename std::tuple_element<Index,StorageType>::type& get()
    {return std::get<Index>(_values);}

    template <Size Index>
    const typename std::tuple_element<Index,StorageType>::type& get() const
    {return std::get<Index>(_values);}

    StorageType& values()
    {return _values;}

    const StorageType& values() const
    {return _values;}

private:
    StorageType _values;
};


} // namespace cie::mp


#endif