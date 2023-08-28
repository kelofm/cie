#ifndef CIE_UTILS_STACK_HPP
#define CIE_UTILS_STACK_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Incldues ---
#include <vector>


namespace cie::utils {


template <class ValueType>
class Stack
{
public:
    using value_type     = ValueType;
    using container_type = std::vector<value_type>;

public:
    Stack( const value_type& r_initializer, Size capacity );

    void push_back( const value_type& r_vector );

    value_type& back();

    value_type& operator[]( size_t index );

private:
    Size           _capacity;
    Size           _counter;
    container_type _values;
};


} // namespace cie::utils

#include "packages/stack/impl/Stack_impl.hpp"

#endif