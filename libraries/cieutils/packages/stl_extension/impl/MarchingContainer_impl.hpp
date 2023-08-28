#ifndef CIE_UTILS_MARCHING_CONTAINER_IMPL_HPP
#define CIE_UTILS_MARCHING_CONTAINER_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::utils {


template <class ValueType>
template <class ContainerType>
requires concepts::Container<ContainerType,ValueType>
MarchingContainer<ValueType>::MarchingContainer( const ContainerType& r_initializer ) :
    _capacity( r_initializer.size() ),
    _counter( r_initializer.size() ),
    _container( r_initializer.size() )
{
    std::copy( r_initializer.begin(),
               r_initializer.end(),
               this->_container.begin() );
}


template <class ValueType>
MarchingContainer<ValueType>::MarchingContainer( Size capacity, const ValueType& r_initializer ) :
    _capacity( capacity ),
    _counter( capacity ),
    _container( capacity, r_initializer )
{
}


template <class ValueType>
MarchingContainer<ValueType>::MarchingContainer( Size capacity ) :
    _capacity( capacity ),
    _counter( 0 ),
    _container( capacity )
{
}


template <class ValueType>
template <class ...Args>
inline void
MarchingContainer<ValueType>::emplace_back( Args&&... args )
{
    this->_counter++;
    this->_container[ this->_counter % this->_capacity ] = std::move(
        ValueType( std::forward<Args>(args)... )
    );
}


template <class ValueType>
inline void
MarchingContainer<ValueType>::push_back( const ValueType& r_value )
{
    this->_counter++;
    this->_container[ this->_counter % this->_capacity ] = r_value;
}


template <class ValueType>
inline ValueType&
MarchingContainer<ValueType>::back()
{
    return this->_container[ this->_counter % this->_capacity ];
}


template <class ValueType>
inline const ValueType&
MarchingContainer<ValueType>::back() const
{
    return this->_container[ this->_counter % this->_capacity ];
}


template <class ValueType>
inline ValueType&
MarchingContainer<ValueType>::operator[]( Size index )
{
    CIE_OUT_OF_RANGE_CHECK( this->outOfRangeCheck(index) )
    return this->_container[ (index+1) % this->_capacity ];
}


template <class ValueType>
inline const ValueType&
MarchingContainer<ValueType>::operator[]( Size index ) const
{
    CIE_OUT_OF_RANGE_CHECK( this->outOfRangeCheck(index) )
    return this->_container[ (index+1) % this->_capacity ];
}


template <class ValueType>
inline Size
MarchingContainer<ValueType>::size() const
{
    return this->_counter;
}


template <class ValueType>
inline bool
MarchingContainer<ValueType>::outOfRangeCheck( Size index )
{
    if ( this->_counter >= this->_capacity ) [[likely]]
        return this->_counter - this->_capacity <= index
               &&
               index < this->_counter;
    else
        return index < this->_counter;
}


} // namespace cie::utils


#endif