#ifndef CIE_UTILITIES_CACHE_IMPL_HPP
#define CIE_UTILITIES_CACHE_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/assertions.hpp"


namespace cie::utils {


template <  class InputType,
            class StoredType >
inline typename AbsCache<InputType,StoredType>::internal_iterator
AbsCache<InputType,StoredType>::insert( const InputType& input,
                                        typename AbsCache::generator_function generator,
                                        bool force )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Hash input and check if it's inserted
    Size id     = this->hash(input);
    auto mapIt  = _map.find(id);

    if (mapIt == _map.end())
        mapIt = _map.emplace( id, generator(input) ).first;
    else if ( force )
        mapIt->second = generator(input);

    return mapIt;

    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline typename AbsCache<InputType,StoredType>::internal_iterator
AbsCache<InputType,StoredType>::insert( const InputType& input,
                                        const StoredType& value,
                                        bool force )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Hash input and check if it's inserted
    Size id     = this->hash(input);
    auto mapIt  = _map.find(id);

    if (mapIt == _map.end())
        mapIt = _map.emplace( id, value ).first;
    else if ( force )
        mapIt->second = value;

    return mapIt;

    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline const typename AbsCache<InputType,StoredType>::stored_type&
AbsCache<InputType,StoredType>::operator[]( Size inputID ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto mapIt = _map.find(inputID);
    if (mapIt == _map.end())
    {
        CIE_THROW(
            OutOfRangeException,
            "Input with ID " + std::to_string(inputID) + " does not have a recorded value in the cache"
        )
    }
    return mapIt->second;

    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline const typename AbsCache<InputType,StoredType>::stored_type&
AbsCache<InputType,StoredType>::operator[]( const InputType& input ) const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return this->operator[]( this->hash(input) );
    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline bool
AbsCache<InputType,StoredType>::cached( Size id ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (_map.find(id) != _map.end())
        return true;
    return false;

    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline bool
AbsCache<InputType,StoredType>::cached( const InputType& input ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return this->cached( this->hash(input) );

    CIE_END_EXCEPTION_TRACING
}


template <  class InputType,
            class StoredType >
inline void
AbsCache<InputType,StoredType>::erase( Size id )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto mapIt = _map.find(id);
    if (mapIt != _map.end())
        _map.erase(mapIt);

    CIE_END_EXCEPTION_TRACING
}


template <  concepts::Container InputContainer,
            concepts::Container OutputContainer >
inline Size
ContainerCache<InputContainer,OutputContainer>::hash( const InputContainer& container ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    Size seed = container.size();
    for ( auto it=container.begin(); it!=container.end(); it++ )
        seed ^= Size(&(*it)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils


#endif