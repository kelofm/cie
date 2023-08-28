#ifndef CIE_LINALG_MATRIX_OPERATORS_HPP
#define CIE_LINALG_MATRIX_OPERATORS_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- Internal Includes ---
#include "packages/matrix/inc/StaticEigenMatrix.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/matrix/inc/concepts.hpp"
#include <packages/matrix/inc/EigenArray.hpp>


namespace cie::linalg {


/* --- Matrix-Scalar Operators --- */


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator+=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs);


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator-=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs);


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
void
operator*=(EigenMatrix<MatrixType>& r_lhs, ValueType rhs);


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
EigenMatrix<MatrixType>
operator*(ValueType lhs, const EigenMatrix<MatrixType>& r_rhs);


template <class MatrixType, class ValueType>
requires std::is_same_v<ValueType, typename MatrixType::value_type>
EigenMatrix<MatrixType>
operator*(const EigenMatrix<MatrixType>& r_lhs, ValueType rhs);


/* --- Matrix-Matrix Operators --- */


template <class MatrixType>
EigenMatrix<MatrixType>
operator+(const EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs);


template <class MatrixType>
EigenMatrix<MatrixType>
operator-(const EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs);


template <class MatrixType>
void
operator+=(EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs);


template <class MatrixType>
void
operator-=(EigenMatrix<MatrixType>& r_lhs, const EigenMatrix<MatrixType>& r_rhs);


template <class TValue>
EigenVector<TValue>
operator*(Ref<const DynamicEigenMatrix<TValue>> r_lhs, Ref<const EigenVector<TValue>> r_rhs);


template <class ValueType>
DynamicEigenMatrix<ValueType>
operator*(const DynamicEigenMatrix<ValueType>& r_lhs, const DynamicEigenMatrix<ValueType>& r_rhs);


template <class ValueType, int LHSSize, int InnerSize, int RHSSize>
StaticEigenMatrix<ValueType, LHSSize, RHSSize>
operator*(const StaticEigenMatrix<ValueType, LHSSize, InnerSize>& r_lhs,
          const StaticEigenMatrix<ValueType, InnerSize, RHSSize>& r_rhs) noexcept;


} // namespace cie::linalg

#include "packages/overloads/impl/matrix_operators_impl.hpp"

#endif
