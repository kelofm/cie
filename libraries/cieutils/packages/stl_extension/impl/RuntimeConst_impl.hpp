#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/RuntimeConst.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::utils {


template <class T>
RuntimeConst<T>::RuntimeConst()
    : _p_value( nullptr ),
      _p_constValue( nullptr )
{
}


template <class T>
RuntimeConst<T>::RuntimeConst(value_type* p_value)
    : _p_value( p_value ),
      _p_constValue( p_value )
{
}


template <class T>
RuntimeConst<T>::RuntimeConst(const value_type* p_value)
    : _p_value( nullptr ),
      _p_constValue( p_value )
{
}


template <class T>
RuntimeConst<T>&
RuntimeConst<T>::operator=(value_type* p_value) {
    _p_value      = p_value;
    _p_constValue = p_value;
    return *this;
}


template <class T>
RuntimeConst<T>&
RuntimeConst<T>::operator=( const value_type* p_value ) const {
    _p_value      = nullptr;
    _p_constValue = p_value;
    return *this;
}


template <class T>
typename RuntimeConst<T>::value_type&
RuntimeConst<T>::operator*() {
    CIE_CHECK_POINTER(_p_value)
    return *this->get();
}


template <class T>
const typename RuntimeConst<T>::value_type&
RuntimeConst<T>::operator*() const {
    CIE_CHECK_POINTER(_p_constValue)
    return *this->get();
}


template <class T>
typename RuntimeConst<T>::value_type*
RuntimeConst<T>::operator->() {
    CIE_CHECK_POINTER(_p_value)
    return this->get();
}


template <class T>
const typename RuntimeConst<T>::value_type*
RuntimeConst<T>::operator->() const {
    CIE_CHECK_POINTER(_p_constValue)
    return this->get();
}


template <class T>
typename RuntimeConst<T>::value_type*
RuntimeConst<T>::get() {
    return _p_value;
}


template <class T>
const typename RuntimeConst<T>::value_type*
RuntimeConst<T>::get() const {
    return _p_constValue;
}


template <class T>
const typename RuntimeConst<T>::value_type*
RuntimeConst<T>::getConst() const {
    return static_cast<const RuntimeConst<T>*>(this)->get();
}


template <class T>
bool
RuntimeConst<T>::isConst() const {
    return _p_constValue != nullptr && _p_value == nullptr;
}


template <class T>
void
RuntimeConst<T>::set(value_type* p_value ) {
    _p_value      = p_value;
    _p_constValue = p_value;
}


template <class T>
void
RuntimeConst<T>::set(const value_type* p_value) {
    _p_value      = nullptr;
    _p_constValue = p_value;
}


template <class T>
RuntimeConst<T>::operator T*() {
    return this->get();
}


template <class T>
RuntimeConst<T>::operator const T*() const {
    return this->get();
}


} // namespace cie::utils
