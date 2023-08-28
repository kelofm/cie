#ifndef CIE_UTILS_UNSTRUCTURED_VIEW_IMPL_HPP
#define CIE_UTILS_UNSTRUCTURED_VIEW_IMPL_HPP

// --- Utility Includes ---
#include "packages/ranges/inc/UnstructuredView.hpp"
#include "packages/types/inc/modifiers.hpp"
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::utils {


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
UnstructuredView<TValue,TContainer>::UnstructuredView(TContainer&& r_pointers)
    : _pointers(std::move(r_pointers))
{
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
template <class TIterator>
requires (concepts::Iterator<TIterator,TValue*> || concepts::Iterator<TIterator,const TValue*>)
UnstructuredView<TValue,TContainer>::UnstructuredView(TIterator begin, TIterator end)
{
    CIE_BEGIN_EXCEPTION_TRACING

    utils::reserve(_pointers, std::distance(begin, end));
    std::copy(begin, end, std::back_inserter(_pointers));

    CIE_END_EXCEPTION_TRACING
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
template <class TIterator>
requires (concepts::Iterator<TIterator,const TValue> || concepts::Iterator<TIterator,TValue>)
UnstructuredView<TValue,TContainer>::UnstructuredView(TIterator begin, TIterator end)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using Value = typename CopyConstQualifier<typename TContainer::value_type,TValue>::Type;

    utils::reserve(_pointers, std::distance(begin, end));
    std::transform(begin, end, std::back_inserter(_pointers), [](Value& ref) {return &ref;});

    CIE_END_EXCEPTION_TRACING
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
typename UnstructuredView<TValue,TContainer>::const_reference
UnstructuredView<TValue,TContainer>::at(size_type index) const
{
    CIE_OUT_OF_RANGE_CHECK(index < this->size())
    return *_pointers.at(index);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
typename UnstructuredView<TValue,TContainer>::reference
UnstructuredView<TValue,TContainer>::at(size_type index)
requires isMutableView
{
    CIE_OUT_OF_RANGE_CHECK(index < this->size())
    return *_pointers.at(index);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
typename UnstructuredView<TValue,TContainer>::const_reference
UnstructuredView<TValue,TContainer>::operator[](size_type index) const
{
    CIE_OUT_OF_RANGE_CHECK(index < this->size())
    return *_pointers.at(index);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
typename UnstructuredView<TValue,TContainer>::reference
UnstructuredView<TValue,TContainer>::operator[](size_type index)
requires isMutableView
{
    CIE_OUT_OF_RANGE_CHECK(index < this->size())
    return *_pointers.at(index);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline void
UnstructuredView<TValue,TContainer>::push_back(ContainerValue& r_value)
{
    _pointers.push_back(&r_value);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline void
UnstructuredView<TValue,TContainer>::emplace_back(ContainerValue& r_value)
{
    _pointers.push_back(&r_value);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline void
UnstructuredView<TValue,TContainer>::reserve(size_type capacity)
{
    _pointers.reserve(capacity);
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline Size
UnstructuredView<TValue,TContainer>::size() const noexcept
{
    return _pointers.size();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline Size
UnstructuredView<TValue,TContainer>::capacity() const noexcept
{
    return _pointers.capacity();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::const_iterator
UnstructuredView<TValue,TContainer>::begin() const noexcept
{
    return _pointers.begin();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::iterator
UnstructuredView<TValue,TContainer>::begin() noexcept
requires isMutableView
{
    return _pointers.begin();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::const_iterator
UnstructuredView<TValue,TContainer>::end() const noexcept
{
    return _pointers.end();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::iterator
UnstructuredView<TValue,TContainer>::end() noexcept
requires isMutableView
{
    return _pointers.end();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::const_reverse_iterator
UnstructuredView<TValue,TContainer>::rbegin() const noexcept
{
    return _pointers.rbegin();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::reverse_iterator
UnstructuredView<TValue,TContainer>::rbegin() noexcept
requires isMutableView
{
    return _pointers.rbegin();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::const_reverse_iterator
UnstructuredView<TValue,TContainer>::rend() const noexcept
{
    return _pointers.rend();
}


template <class TValue, class TContainer>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
inline typename UnstructuredView<TValue,TContainer>::reverse_iterator
UnstructuredView<TValue,TContainer>::rend() noexcept
requires isMutableView
{
    return _pointers.rend();
}


} // namespace cie::utils

#endif
