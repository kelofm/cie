#ifndef CIE_UTILS_STL_EXTENSION_RUNTIME_CONST_HPP
#define CIE_UTILS_STL_EXTENSION_RUNTIME_CONST_HPP


namespace cie::utils {


template <class ValueType>
class RuntimeConst
{
public:
    using value_type = ValueType;

public:
    RuntimeConst();

    RuntimeConst( value_type* p_value );

    RuntimeConst( const value_type* p_value );

    RuntimeConst<ValueType>& operator=( value_type* p_value );

    RuntimeConst<ValueType>& operator=( const value_type* p_value ) const;

    value_type& operator*();

    const value_type& operator*() const;

    value_type* operator->();

    const value_type* operator->() const;

    value_type* get();

    const value_type* get() const;

    const value_type* getConst() const;

    void set( value_type* p_value );

    void set( const value_type* p_value );

    bool isConst() const;

    explicit operator ValueType*();

    operator const ValueType*() const;

private:
    value_type*       _p_value;
    const value_type* _p_constValue;
}; // class RuntimeConst


} // namespace cie::utils

#include "packages/stl_extension/impl/RuntimeConst_impl.hpp"

#endif