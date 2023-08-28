#ifndef CIE_LINALG_MATRIX_OPERATORS_IMPL_HPP
#define CIE_LINALG_MATRIX_OPERATORS_IMPL_HPP

// --- LinAlg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include <packages/matrix/inc/DynamicEigenMatrix.hpp>


namespace cie::linalg {


/* --- Matrix-Scalar Operators --- */


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator+=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_lhs.wrapped() += rhs;

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator-=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_lhs.wrapped() -= rhs;

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator*=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_lhs.wrapped() *= rhs;

    CIE_END_EXCEPTION_TRACING
}


/* --- Matrix-Matrix Operators --- */


template <class MatrixType>
void
operator+=(EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_lhs.wrapped() += r_rhs.wrapped();

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType>
void
operator-=(EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_lhs.wrapped() -= r_rhs.wrapped();

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType>
EigenMatrix<MatrixType>
operator+(const EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return EigenMatrix<MatrixType>((r_lhs.wrapped() + r_rhs.wrapped()).eval());

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType>
EigenMatrix<MatrixType>
operator-(const EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return EigenMatrix<MatrixType>((r_lhs.wrapped() - r_rhs.wrapped()).eval());

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
EigenMatrix<MatrixType>
operator*(ValueType lhs, const EigenMatrix<MatrixType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return EigenMatrix<MatrixType>(
        lhs * r_rhs.wrapped()
    );

    CIE_END_EXCEPTION_TRACING
}


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
EigenMatrix<MatrixType>
operator*(const EigenMatrix<MatrixType>& r_lhs, ValueType rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return EigenMatrix<MatrixType>(
        r_lhs.wrapped() * rhs
    );

    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
EigenVector<TValue>
operator*(Ref<const DynamicEigenMatrix<TValue>> r_lhs, Ref<const EigenVector<TValue>> r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return EigenVector<TValue>(
        (r_lhs.wrapped() * r_rhs.wrapped()).eval()
    );
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType>
DynamicEigenMatrix<ValueType>
operator*(const DynamicEigenMatrix<ValueType>& r_lhs, const DynamicEigenMatrix<ValueType>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return DynamicEigenMatrix<ValueType>(
        (r_lhs.wrapped() * r_rhs.wrapped()).eval()
    );
    CIE_END_EXCEPTION_TRACING
}


template <class ValueType, int LHSSize, int InnerSize, int RHSSize>
StaticEigenMatrix<ValueType, LHSSize, RHSSize>
operator*(const StaticEigenMatrix<ValueType, LHSSize, InnerSize>& r_lhs,
           const StaticEigenMatrix<ValueType, InnerSize, RHSSize>& r_rhs) noexcept
{
    return StaticEigenMatrix<ValueType, LHSSize, RHSSize>(
        (r_lhs.wrapped() * r_rhs.wrapped())
    );
}


} // namespace cie::linalg


#endif