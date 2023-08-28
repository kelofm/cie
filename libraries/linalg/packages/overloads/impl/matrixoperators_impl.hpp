#ifndef LINALG_MATRIX_OPERATORS_IMPL_HPP
#define LINALG_MATRIX_OPERATORS_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STD Incldues ---
#include <string>
#include <algorithm>

namespace cie::linalg {


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator+(const Matrix<ValueType>& matrix, ScalarType scalar)
{
    std::vector<ValueType> rowMajorData(matrix.size1() * matrix.size2());
    std::transform(
        matrix.data().begin(),
        matrix.data().end(),
        rowMajorData.begin(),
        [scalar](auto matrixIt) -> ScalarType
            {return matrixIt + scalar;}
        );

    return Matrix<ValueType>(rowMajorData, matrix.size1());
}


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator+(ScalarType scalar, const Matrix<ValueType>& matrix)
{
    return matrix + scalar;
}


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator-(const Matrix<ValueType>& matrix, ScalarType scalar)
{
    return matrix + (-scalar);
}


template <concepts::Numeric ValueType>
Matrix<ValueType> operator+(const Matrix<ValueType>& lhs, const Matrix<ValueType>& rhs)
{
    // Check sizes
    if (lhs.size1() != rhs.size1()  ||  lhs.size2() != rhs.size2())
        CIE_THROW( OutOfRangeException, "Matrix size mismatch!" )

    std::vector<ValueType> rowMajorData(lhs.size1()*lhs.size2());
    std::transform(
        lhs.data().begin(),
        lhs.data().end(),
        rhs.data().begin(),
        rowMajorData.begin(),
        [](auto lhsIt, auto rhsIt) -> ValueType
            {return lhsIt + rhsIt;}
    );

    return Matrix<ValueType>(rowMajorData,lhs.size1());
}


template <concepts::Numeric ValueType>
Matrix<ValueType> operator-(const Matrix<ValueType>& lhs, const Matrix<ValueType>& rhs)
{
    // Check sizes
    if (lhs.size1() != rhs.size1()  ||  lhs.size2() != rhs.size2())
        CIE_THROW( OutOfRangeException, "Matrix size mismatch!" )

    std::vector<ValueType> rowMajorData(lhs.size1()*lhs.size2());
    std::transform(
        lhs.data().begin(),
        lhs.data().end(),
        rhs.data().begin(),
        rowMajorData.begin(),
        [](auto lhsIt, auto rhsIt) -> ValueType
            {return lhsIt - rhsIt;}
    );

    return Matrix<ValueType>(rowMajorData,lhs.size1());
}


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator*(const Matrix<ValueType>& matrix, ScalarType scalar)
{
    std::vector<ValueType> rowMajorData(matrix.size1() * matrix.size2());
    std::transform(
        matrix.data().begin(),
        matrix.data().end(),
        rowMajorData.begin(),
        [scalar](auto matrixIt) -> ScalarType
            {return matrixIt * scalar;}
        );

    return Matrix<ValueType>(rowMajorData, matrix.size1());
}


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator*(ScalarType scalar, const Matrix<ValueType>& matrix)
{
    return matrix * scalar;
}


template <concepts::Numeric ValueType, concepts::Numeric ScalarType>
Matrix<ValueType> operator/(const Matrix<ValueType>& matrix, ScalarType scalar)
{
    if (std::abs(scalar)==0)
        CIE_THROW( DivisionByZeroException, "" )
    return matrix * (1/scalar);
}


template <concepts::Numeric ValueType, concepts::NumericContainer ContainerType>
ContainerType operator*(const ContainerType& vector, const Matrix<ValueType>& matrix)
{
    // Check size compatibility
    if (vector.size() != matrix.size1())
        CIE_THROW(
            OutOfRangeException,
            "Cannot multiply vector of size " +
            std::to_string(vector.size()) +
            " by matrix of size " +
            std::to_string(matrix.size1()) + "x" + std::to_string(matrix.size2())
        )

    ContainerType result(matrix.size2());
    std::fill( result.begin(), result.end(), 0 );

    for (Size j=0; j<matrix.size2(); ++j){
        for (Size i=0; i<matrix.size1(); ++i)
            result[j] += vector[i] * matrix(i,j);
    }

    return result;
}


template <concepts::Numeric ValueType, Size N>
StaticArray<ValueType,N> operator*(const StaticArray<ValueType,N>& vector, const Matrix<ValueType>& matrix)
{
    // Check size compatibility
    if (vector.size() != matrix.size1())
        CIE_THROW(
            OutOfRangeException,
            "Cannot multiply vector of size " +
            std::to_string(vector.size()) +
            " by matrix of size " +
            std::to_string(matrix.size1()) + "x" + std::to_string(matrix.size2())
        )

    StaticArray<ValueType,N> result;
    std::fill( result.begin(), result.end(), 0 );

    for (Size j=0; j<matrix.size2(); ++j){
        for (Size i=0; i<matrix.size1(); ++i)
            result[j] += vector[i] * matrix(i,j);
    }

    return result;
}


template <concepts::Numeric ValueType, concepts::NumericContainer ContainerType>
ContainerType operator*(const Matrix<ValueType>& matrix, const ContainerType& vector)
{
    // Check size compatibility
    if (vector.size() != matrix.size2())
        CIE_THROW(
            OutOfRangeException,
            "Cannot multiply matrix of size " +
            std::to_string(matrix.size1()) + "x" + std::to_string(matrix.size2()) +
            " by vector of size " +
            std::to_string(vector.size())
        )

    ContainerType result(matrix.size1());
    std::fill( result.begin(), result.end(), 0 );

    for (Size i=0; i<matrix.size1(); ++i){
        for (Size j=0; j<matrix.size2(); ++j)
            result[i] += matrix(i,j) * vector[j];
    }

    return result;
}


template <concepts::Numeric ValueType, Size N>
StaticArray<ValueType,N> operator*(const Matrix<ValueType>& matrix, const StaticArray<ValueType,N>& vector)
{
    // Check size compatibility
    if (vector.size() != matrix.size2())
        CIE_THROW(
            OutOfRangeException,
            "Cannot multiply matrix of size " +
            std::to_string(matrix.size1()) + "x" + std::to_string(matrix.size2()) +
            " by vector of size " +
            std::to_string(vector.size())
        )

    StaticArray<ValueType,N> result;
    std::fill( result.begin(), result.end(), 0 );

    for (Size i=0; i<matrix.size1(); ++i){
        for (Size j=0; j<matrix.size2(); ++j)
            result[i] += matrix(i,j) * vector[j];
    }

    return result;
}


template <concepts::Numeric ValueType>
Matrix<ValueType> operator*(const Matrix<ValueType>& lhs, const Matrix<ValueType>& rhs)
{
    // Check size compatibility
    if (lhs.size2() != rhs.size1())
        CIE_THROW(
            OutOfRangeException,
            "Cannot multiply matrix of size " +
            std::to_string(lhs.size1()) + "x" + std::to_string(lhs.size2()) +
            " by matrix of size " +
            std::to_string(rhs.size1()) + "x" + std::to_string(rhs.size2())
        )

    Matrix<ValueType> result(lhs.size1(), rhs.size2(), 0.0);

    for (Size i=0; i<lhs.size1(); ++i){
        for (Size j=0; j<rhs.size2(); ++j){

            for (Size k=0; k<lhs.size2(); ++k){
                result(i,j) += lhs(i,k) * rhs(k,j);
            }

        }
    }

    return result;
}


} // namespace cie::linalg

#endif