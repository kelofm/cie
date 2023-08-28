#ifndef CIE_UTILS_STL_EXTENSION_STATE_ITERATOR_IMPL_HPP
#define CIE_UTILS_STL_EXTENSION_STATE_ITERATOR_IMPL_HPP

// --- Internal Includes ---
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>
#include <functional>
#include <vector>


namespace cie::utils {


template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
StateIterator<IteratorType>::StateIterator( IteratorType begin,
                                            IteratorType end ) :
    _begin(begin),
    _end(end)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->reset();

    CIE_END_EXCEPTION_TRACING
}


template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
inline void
StateIterator<IteratorType>::reset()
{
    CIE_BEGIN_EXCEPTION_TRACING

    utils::resize( _state, std::distance(_begin,_end) );
    std::transform( _begin,
                    _end,
                    _state.begin(),
                    [](const auto& container) -> subiterator_type
                    { return container.begin(); } );

    CIE_END_EXCEPTION_TRACING
}


template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
inline StateIterator<IteratorType>&
StateIterator<IteratorType>::operator++()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK( _state.size() == Size(std::distance(_begin,_end)) )

    auto stateIt        = _state.begin();
    auto subContainerIt = _begin;

    while( subContainerIt != _end )
    {
        if ( (*stateIt)+1 != subContainerIt->end() )    // Found the place to increment
        {
            ++(*stateIt);
            break;
        }
        else                                            // Overflow
        {
            *stateIt = subContainerIt->begin();
            stateIt++;
            subContainerIt++;
        }
    } // while incrementable

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
inline const typename StateIterator<IteratorType>::state_container&
StateIterator<IteratorType>::operator*() const
{
    return _state;
}


template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
StateIterator<IteratorType>::StateIterator()
{
}




template <concepts::Container ContainerType>
InternalStateIterator<ContainerType>::InternalStateIterator( const ContainerType& container,
                                                             Size size ) :
    StateIterator<typename std::vector<detail::container_reference_wrapper<const ContainerType>>::const_iterator>()
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (Size i=0; i<size; ++i)
        _dummyContainer.emplace_back(container);

    this->_begin    = _dummyContainer.begin();
    this->_end      = _dummyContainer.end();
    this->reset();

    CIE_END_EXCEPTION_TRACING
}




// ---------------------------------------------------------
// CONVENIENCE FUNCTIONS
// ---------------------------------------------------------

template <class IteratorType>
StateIterator<IteratorType>
makeStateIterator(  IteratorType begin,
                    IteratorType end )
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
{
    CIE_BEGIN_EXCEPTION_TRACING

    return StateIterator<IteratorType>(begin,end);

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container ContainerType>
StateIterator<typename ContainerType::const_iterator>
makeStateIterator( const ContainerType& container )
requires concepts::Container<typename ContainerType::value_type>
{
    CIE_BEGIN_EXCEPTION_TRACING

    return StateIterator<typename ContainerType::const_iterator>( container.begin(), container.end() );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container ContainerType>
InternalStateIterator<const ContainerType>
makeInternalStateIterator(const ContainerType& container,
                          Size size )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return InternalStateIterator<const ContainerType>( container, size );

    CIE_END_EXCEPTION_TRACING
}





} // namespace cie::utils


#endif