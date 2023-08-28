#ifndef CIE_UTILS_STACK_IMPL_HPP
#define CIE_UTILS_STACK_IMPL_HPP

// --- STL Includes ---
#include <algorithm>


namespace cie::utils {


template <class ValueType>
Stack<ValueType>::Stack( const ValueType& r_initializer, Size capacity ) :
    _capacity( capacity ),
    _counter( capacity ),
    _values( capacity )
{
    std::fill(
        this->_values.begin(),
        this->_values.end(),
        r_initializer
    );
}


template <class ValueType>
inline void
Stack<ValueType>::push_back( const ValueType& r_value )
{
    (this->_counter)++;
    this->_values[ this->_counter % this->_capacity ] = r_value;
}


template <class ValueType>
inline ValueType&
Stack<ValueType>::back()
{
    return this->_values[ this->_counter % this->_capacity ];
}


template <class ValueType>
inline ValueType&
Stack<ValueType>::operator[]( Size index )
{
    return this->_values[ (index+1) % this->_capacity ];
}


} // namespace cie::utils


#endif