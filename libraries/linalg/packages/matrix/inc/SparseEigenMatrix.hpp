#ifndef CIE_LINALG_SPARSE_EIGEN_MATRIX_HPP
#define CIE_LINALG_SPARSE_EIGEN_MATRIX_HPP

// --- External Includes ---
#include <Eigen/Sparse>

// --- Internal Includes ---
#include "packages/matrix/inc/EigenMatrix.hpp"


namespace cie::linalg {


template <class ValueType>
using SparseEigenMatrix = EigenMatrix<
    Eigen::SparseMatrix<ValueType>
>;


} // namespace cie::linalg


#endif