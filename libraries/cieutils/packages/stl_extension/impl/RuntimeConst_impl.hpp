#ifndef CIE_UTILS_STL_EXTENSION_RUNTIME_CONST_IMPL_HPP
#define CIE_UTILS_STL_EXTENSION_RUNTIME_CONST_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"


namespace cie::utils {


template <class ValueType>
RuntimeConst<ValueType>::RuntimeConst()
    : _p_value( nullptr ),
      _p_constValue( nullptr )
{
}


template <class ValueType>
RuntimeConst<ValueType>::RuntimeConst( typename RuntimeConst<ValueType>::value_type* p_value )
    : _p_value( p_value ),
      _p_constValue( p_value )
{
}


template <class ValueType>
RuntimeConst<ValueType>::RuntimeConst( const RuntimeConst<ValueType>::value_type* p_value )
    : _p_value( nullptr ),
      _p_constValue( p_value )
{
}


template <class ValueType>
inline RuntimeConst<ValueType>&
RuntimeConst<ValueType>::operator=( typename RuntimeConst<ValueType>::value_type* p_value )
{
    _p_value      = p_value;
    _p_constValue = p_value;

    return *this;
}


template <class ValueType>
inline RuntimeConst<ValueType>&
RuntimeConst<ValueType>::operator=( const typename RuntimeConst<ValueType>::value_type* p_value ) const
{
    _p_value      = nullptr;
    _p_constValue = p_value;

    return *this;
}


template <class ValueType>
inline typename RuntimeConst<ValueType>::value_type&
RuntimeConst<ValueType>::operator*()
{
    CIE_BEGIN_EXCEPTION_TRACING
    return *this->get();
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
inline const typename RuntimeConst<ValueType>::value_type&
RuntimeConst<ValueType>::operator*() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return *this->get();
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
inline typename RuntimeConst<ValueType>::value_type*
RuntimeConst<ValueType>::operator->()
{
    CIE_BEGIN_EXCEPTION_TRACING
    return this->get();
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
inline const typename RuntimeConst<ValueType>::value_type*
RuntimeConst<ValueType>::operator->() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return this->get();
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
inline typename RuntimeConst<ValueType>::value_type*
RuntimeConst<ValueType>::get()
{
    CIE_CHECK(
        _p_value != nullptr,
        "attempt to get non-const pointer to const object"
    )

    return _p_value;
}


template <class ValueType>
inline const typename RuntimeConst<ValueType>::value_type*
RuntimeConst<ValueType>::get() const
{
    CIE_CHECK_POINTER( _p_constValue )

    return _p_constValue;
}


template <class ValueType>
inline const typename RuntimeConst<ValueType>::value_type*
RuntimeConst<ValueType>::getConst() const
{
    return const_cast<const RuntimeConst<ValueType>*>(this)->get();
}


template <class ValueType>
inline bool
RuntimeConst<ValueType>::isConst() const
{
    CIE_CHECK_POINTER( _p_constValue )
    return _p_value == nullptr;
}


template <class ValueType>
inline void
RuntimeConst<ValueType>::set( typename RuntimeConst<ValueType>::value_type* p_value )
{
    _p_value      = p_value;
    _p_constValue = p_value;
}


template <class ValueType>
inline void
RuntimeConst<ValueType>::set( const typename RuntimeConst<ValueType>::value_type* p_value )
{
    _p_value      = nullptr;
    _p_constValue = p_value;
}


template <class ValueType>
RuntimeConst<ValueType>::operator ValueType*()
{
    CIE_BEGIN_EXCEPTION_TRACING
    return this->get();
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
RuntimeConst<ValueType>::operator const ValueType*() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return this->get();
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils

#endif