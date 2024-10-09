#ifndef CIE_UTILS_THREAD_POOL_STORAGE_HPP
#define CIE_UTILS_THREAD_POOL_STORAGE_HPP

// --- External Includes ---
#include "tsl/robin_map.h"

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadStorage.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <iterator>
#include <utility>


namespace cie::mp {


/// @ingroup cieutils
template <class TStorage>
class ThreadPoolStorage
{
private:
    using StorageContainer = tsl::robin_map<std::thread::id,TStorage>;

    template <class TValue>
    class Iterator
    {
    public:
        using value_type = TValue;

        using reference = Ref<TValue>;

        using const_reference = Ref<const TValue>;

        using pointer = Ptr<TValue>;

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::forward_iterator_tag;

    private:
        template <class T>
        struct SelectIterator
        {using Type = typename StorageContainer::iterator;};

        template <concepts::Const T>
        struct SelectIterator<T>
        {using Type = typename StorageContainer::const_iterator;};

        using Wrapped = typename SelectIterator<TValue>::Type;

    public:
        Iterator() = delete;

        Iterator(Wrapped wrapped) : _wrapped(wrapped) {}

        Iterator(RightRef<Iterator> r_rhs) noexcept = default;

        Iterator(Ref<const Iterator> r_rhs) = default;

        Iterator& operator=(Iterator&& r_rhs) noexcept = default;

        Iterator& operator=(const Iterator& r_rhs) = default;

        const_reference operator*() const
        {return _wrapped->second;}

        reference operator*()
        requires concepts::NonConst<TValue>
        {return _wrapped.value();}

        pointer operator->() const
        requires concepts::Const<TValue>
        {return &_wrapped->second;}

        pointer operator->()
        requires concepts::NonConst<TValue>
        {return _wrapped.value();}

        Ref<Iterator> operator++() {++_wrapped; return *this;}

        Iterator operator++(int) {Iterator copy(*this); ++(*this); return copy;}

        friend difference_type operator-(Iterator lhs, Iterator rhs) {return lhs._wrapped - rhs._wrapped;}

        friend Bool operator<(Iterator lhs, Iterator rhs) {return lhs._wrapped < rhs._wrapped;}

        friend Bool operator!=(Iterator lhs, Iterator rhs) {return lhs._wrapped != rhs._wrapped;}

    private:
        Wrapped _wrapped;
    }; // class Iterator

public:
    using iterator = Iterator<TStorage>;

    using const_iterator = Iterator<const TStorage>;

    using value_type = TStorage;

    using size_type = typename StorageContainer::size_type;

public:
    ThreadPoolStorage(Ref<ThreadPoolBase> r_pool);

    ThreadPoolStorage(Ref<ThreadPoolBase> p_pool, Ref<const TStorage> r_storage);

    ThreadPoolStorage(ThreadPoolStorage&& r_rhs) noexcept = default;

    ThreadPoolStorage(const ThreadPoolStorage& r_rhs) = default;

    ThreadPoolStorage& operator=(ThreadPoolStorage&& r_rhs) noexcept = default;

    ThreadPoolStorage& operator=(const ThreadPoolStorage& r_rhs) = default;

    Ref<TStorage> getLocalStorage();

    Ref<const TStorage> getLocalStorage() const;

    void setLocalStorage(RightRef<TStorage> r_storage);

    void setLocalStorage(Ref<const TStorage> r_storage);

    Ref<value_type> at(std::thread::id id) noexcept
    {return _storages[id];}

    Ref<const value_type> at(std::thread::id id) const noexcept
    {return _storages[id];}

    Ref<value_type> operator[](std::thread::id id) noexcept
    {return _storages[id];}

    Ref<const value_type> operator[](std::thread::id id) const noexcept
    {return _storages[id];}

    Size size() const noexcept
    {return _storages.size();}

    iterator begin() noexcept
    {return iterator(_storages.begin());}

    const_iterator begin() const noexcept
    {return const_iterator(_storages.begin());}

    iterator end() noexcept
    {return iterator(_storages.end());}

    const_iterator end() const noexcept
    {return const_iterator(_storages.end());}

private:
    std::reference_wrapper<ThreadPoolBase> _r_pool;

    StorageContainer _storages;
}; // class ThreadPoolStorage


} // namespace cie::mp

#include "packages/concurrency/impl/ThreadPoolStorage_impl.hpp"

#endif