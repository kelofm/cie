#ifndef CIE_UTILS_STL_EXTENSION_ZIP_IMPL_HPP
#define CIE_UTILS_STL_EXTENSION_ZIP_IMPL_HPP


#include <utility>
namespace cie::utils {


namespace detail {


template <class ...TContainers>
template <class TValue>
inline TValue
ZipIterator<TContainers...>::makeRefTuple()
{
    return this->makeRefTupleImpl<TValue>(
        std::make_index_sequence<sizeof...(TContainers)>()
    );
}


template <class ...TContainers>
template <class TValue, std::size_t ...Indices>
inline TValue
ZipIterator<TContainers...>::makeRefTupleImpl(std::index_sequence<Indices...>)
{
    return TValue(*std::get<Indices>(_it)...);
}


template <class ...TContainers>
template <std::size_t ...Indices>
inline typename ZipObject<TContainers...>::Iterator
ZipObject<TContainers...>::beginImpl(std::index_sequence<Indices...>)
{
    using Iterator = typename ZipObject<TContainers...>::Iterator;
    return Iterator(
        typename Iterator::IteratorTuple(std::get<Indices>(_containers).begin()...)
    );
}


template <class ...TContainers>
template <std::size_t ...Indices>
inline typename ZipObject<TContainers...>::Iterator
ZipObject<TContainers...>::endImpl(std::index_sequence<Indices...>)
{
    auto it = this->begin();
    it += _size;
    return it;
}


} // namespace detail


} // namespace cie::utils


#endif