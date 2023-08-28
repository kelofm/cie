#ifndef CIE_LINALG_OVERLOADS_VECTOR_OPERATORS_IMPL_HPP
#define CIE_LINALG_OVERLOADS_VECTOR_OPERATORS_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include <packages/macros/inc/checks.hpp>


namespace cie::linalg {


// -----------------------------------------------------------
// IN-PLACE OPERATORS
// -----------------------------------------------------------

template <class EigenType, class BaseType, class ValueType>
requires (std::is_same_v<ValueType, typename BaseType::value_type> && concepts::NonConst<BaseType>)
void
operator*=( AbsEigenVector<EigenType,BaseType>& r_lhs, ValueType rhs )
{
    r_lhs.wrapped() *= rhs;
}


template <class EigenType, class BaseType>
requires concepts::NonConst<BaseType>
void
operator+=( AbsEigenVector<EigenType, BaseType>& r_lhs,
            const AbsEigenVector<EigenType, BaseType>& r_rhs )
{
    #ifdef __GNUG__
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Warray-bounds"
    #endif
    r_lhs.wrapped() += r_rhs.wrapped();
    #ifdef __GNUG__
    #pragma GCC diagnostic pop
    #endif
}


template <class EigenType, class BaseType>
requires concepts::NonConst<BaseType>
void
operator-=( AbsEigenVector<EigenType, BaseType>& r_lhs,
            const AbsEigenVector<EigenType, BaseType>& r_rhs )
{
    #ifdef __GNUG__
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Warray-bounds"
    #endif
    r_lhs.wrapped() -= r_rhs.wrapped();
    #ifdef __GNUG__
    #pragma GCC diagnostic pop
    #endif
}


// -----------------------------------------------------------
// SCALAR-VECTOR OPERATORS
// -----------------------------------------------------------

template <class EigenType, class BaseType, concepts::Numeric ValueType>
requires concepts::DynamicArray<EigenType>
inline EigenVector<ValueType>
operator*( ValueType lhs, const AbsEigenVector<EigenType,BaseType>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    EigenVector<ValueType> output( r_rhs );
    output.wrapped() *= lhs;

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <class EigenType, class BaseType, concepts::Numeric ValueType>
requires concepts::DynamicArray<EigenType>
inline EigenVector<ValueType>
operator*( const AbsEigenVector<EigenType,BaseType>& r_lhs, ValueType rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    EigenVector<ValueType> output( r_lhs );
    output.wrapped() *= rhs;

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <class EigenType, class BaseType, concepts::Numeric ValueType>
requires concepts::StaticArray<EigenType>
inline EigenArray<ValueType, AbsEigenVector<EigenType,BaseType>::RowTag>
operator*( ValueType lhs, const AbsEigenVector<EigenType,BaseType>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    EigenArray<ValueType, AbsEigenVector<EigenType,BaseType>::RowTag> output( r_rhs );
    output.wrapped() *= lhs;

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <class EigenType, class BaseType, concepts::Numeric ValueType>
requires concepts::StaticArray<EigenType>
inline EigenArray<ValueType, AbsEigenVector<EigenType,BaseType>::RowTag>
operator*( const AbsEigenVector<EigenType,BaseType>& r_lhs, ValueType rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    EigenArray<ValueType, AbsEigenVector<EigenType,BaseType>::RowTag> output( r_lhs );
    output.wrapped() *= rhs;

    return output;

    CIE_END_EXCEPTION_TRACING
}


// -----------------------------------------------------------
// VECTOR-VECTOR OPERATORS
// -----------------------------------------------------------

template <class EigenType, class ValueType, Size ArraySize>
EigenArray<ValueType, ArraySize>
operator+( const AbsEigenVector<EigenType, StaticArray<ValueType,ArraySize>>& r_lhs,
           const AbsEigenVector<EigenType, StaticArray<ValueType,ArraySize>>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_lhs.size() == r_rhs.size())
    EigenArray<ValueType,ArraySize> output = r_lhs.base();
    output.wrapped() += r_rhs.wrapped();

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <class EigenType, class ValueType, Size ArraySize>
EigenArray<ValueType, ArraySize>
operator-( const AbsEigenVector<EigenType, StaticArray<ValueType,ArraySize>>& r_lhs,
           const AbsEigenVector<EigenType, StaticArray<ValueType,ArraySize>>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_lhs.size() == r_rhs.size())
    EigenArray<ValueType,ArraySize> output = r_lhs.base();
    output.wrapped() -= r_rhs.wrapped();

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <class EigenType, class ValueType, class ...Arguments>
EigenVector<ValueType,Arguments...>
operator+( const AbsEigenVector<EigenType, std::vector<ValueType,Arguments...>>& r_lhs,
           const AbsEigenVector<EigenType, std::vector<ValueType,Arguments...>>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_lhs.size() == r_rhs.size())
    EigenVector<ValueType,Arguments...> output = r_lhs.base();
    output.wrapped() += r_rhs.wrapped();

    return output;

    CIE_END_EXCEPTION_TRACING
}

template <class EigenType, class ValueType, class ...Arguments>
EigenVector<ValueType,Arguments...>
operator-( const AbsEigenVector<EigenType, std::vector<ValueType,Arguments...>>& r_lhs,
           const AbsEigenVector<EigenType, std::vector<ValueType,Arguments...>>& r_rhs )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(r_lhs.size() == r_rhs.size())
    EigenVector<ValueType,Arguments...> output = r_lhs.base();
    output.wrapped() -= r_rhs.wrapped();

    return output;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::linalg


#endif