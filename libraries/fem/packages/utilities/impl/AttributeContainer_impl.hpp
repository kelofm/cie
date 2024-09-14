#ifndef CIE_FEM_ATTRIBUTE_CONTAINER_IMPL_HPP
#define CIE_FEM_ATTRIBUTE_CONTAINER_IMPL_HPP

// --- FEM Includes ---
#include "packages/utilities/inc/AttributeContainer.hpp"


namespace cie::fem {


template <class ...TValues>
AttributeContainer<TValues...>::AttributeContainer(Size size)
    : _containers(DynamicArray<TValues>(size) ...)
{
}


namespace Impl {


template <class TIndexSequence, class ...TValues>
struct AttributeContainer {};


template <std::size_t ...Indices, class ...TValues>
struct AttributeContainer<std::index_sequence<Indices...>, TValues...>
{
    using Tuple = std::tuple<DynamicArray<TValues>...>;

    using IndexSequence = std::index_sequence<Indices...>;

    static void resize(Ref<Tuple> rContainers, Size size)
    {
        (..., std::get<Indices>(rContainers).resize(size));
    }

    static void erase(Ref<Tuple> rContainers, Size iBegin, Size iEnd)
    {
        (..., std::get<Indices>(rContainers).erase(
                std::get<Indices>(rContainers).begin() + iBegin,
                std::get<Indices>(rContainers).begin() + iEnd
        ));
    }

    static void clear(Ref<Tuple> rContainers)
    {
        (..., std::get<Indices>(rContainers).clear());
    }

    static void reserve(Ref<Tuple> rContainers, Size capacity)
    {
        (..., std::get<Indices>(rContainers).reserve(capacity));
    }

    static void push_back(Ref<Tuple> rContainers, TValues... values)
    {
        (..., std::get<Indices>(rContainers).push_back(values));
    }
}; // struct AttributeContainer

} // namespace Impl


template <class ...TValues>
void
AttributeContainer<TValues...>::resize(Size size)
{
    Impl::AttributeContainer<std::remove_const_t<decltype(TupleIndices)>,TValues...>::resize(_containers, size);
}


template <class ...TValues>
void
AttributeContainer<TValues...>::erase(Size iBegin, Size iEnd)
{
    Impl::AttributeContainer<std::remove_const_t<decltype(TupleIndices)>,TValues...>::erase(_containers, iBegin, iEnd);
}


template <class ...TValues>
void
AttributeContainer<TValues...>::erase(Size index)
{
    this->erase(index, index + 1);
}


template <class ...TValues>
void
AttributeContainer<TValues...>::clear()
{
    Impl::AttributeContainer<std::remove_const_t<decltype(TupleIndices)>,TValues...>::clear(_containers);
}


template <class ...TValues>
void
AttributeContainer<TValues...>::reserve(Size capacity)
{
    Impl::AttributeContainer<std::remove_const_t<decltype(TupleIndices)>,TValues...>::reserve(_containers, capacity);
}


template <class ...TValues>
void
AttributeContainer<TValues...>::push_back(TValues... values)
{
    Impl::AttributeContainer<std::remove_const_t<decltype(TupleIndices)>,TValues...>::push_back(_containers, values...);
}


template <class ...TValues>
Size
AttributeContainer<TValues...>::size() const noexcept
{
    return std::get<0>(_containers).size();
}


template <class ...TValues>
Size
AttributeContainer<TValues...>::capacity() const noexcept
{
    return std::get<0>(_containers).capacity();
}


template <class ...TValues>
bool
AttributeContainer<TValues...>::empty() const noexcept
{
    return std::get<0>(_containers).empty();
}


template <class ...TValues>
template <class TValue>
TValue
AttributeContainer<TValues...>::at(Size index) const noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).at(index);
}


template <class ...TValues>
template <class TValue>
Ref<TValue>
AttributeContainer<TValues...>::at(Size index) noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).at(index);
}


template <class ...TValues>
std::tuple<TValues...>
AttributeContainer<TValues...>::get(Size index) const noexcept
{
    return std::make_tuple(this->at<TValues>(index)...);
}


template <class ...TValues>
template <class TValue>
typename DynamicArray<TValue>::const_iterator
AttributeContainer<TValues...>::begin() const noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).begin();
}


template <class ...TValues>
template <class TValue>
typename DynamicArray<TValue>::iterator
AttributeContainer<TValues...>::begin() noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).begin();
}


template <class ...TValues>
template <class TValue>
typename DynamicArray<TValue>::const_iterator
AttributeContainer<TValues...>::end() const noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).end();
}


template <class ...TValues>
template <class TValue>
typename DynamicArray<TValue>::iterator
AttributeContainer<TValues...>::end() noexcept
{
    return std::get<DynamicArray<TValue>>(_containers).end();
}


} // namespace cie::fem


#endif
