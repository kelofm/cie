#ifndef CIE_LINALG_STATIC_EIGEN_MATRIX_HPP
#define CIE_LINALG_STATIC_EIGEN_MATRIX_HPP

// --- External Includes ---
#include <Eigen/Dense>

// --- Internal Includes ---
#include "packages/matrix/inc/EigenMatrix.hpp"


namespace cie::linalg {


template <class ValueType, int RowSize, int ColumnSize>
requires (0<RowSize && 0<ColumnSize)
using StaticEigenMatrix = EigenMatrix<
    Eigen::Matrix<ValueType,
                  RowSize,
                  ColumnSize>
>;


} // namespace cie::linalg


#endif