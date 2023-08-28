#ifndef CIE_LINALG_DYNAMIC_EIGEN_MATRIX_HPP
#define CIE_LINALG_DYNAMIC_EIGEN_MATRIX_HPP

// --- External Includes ---
#include <Eigen/Dense>

// --- Internal Includes ---
#include "packages/matrix/inc/EigenMatrix.hpp"


namespace cie::linalg {


template <class ValueType>
using DynamicEigenMatrix = EigenMatrix<
    Eigen::Matrix< ValueType,
                   Eigen::Dynamic,
                   Eigen::Dynamic>
>;


} // namespace cie::linalg


#endif