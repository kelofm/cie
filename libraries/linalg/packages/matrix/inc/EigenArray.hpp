#ifndef CIE_LINALG_EIGEN_ARRAY_HPP
#define CIE_LINALG_EIGEN_ARRAY_HPP

// --- External Includes ---
#include "Eigen/Dense"

// --- LinAlg Includes ---
#include "packages/matrix/inc/EigenMatrix.hpp"


namespace cie::linalg {


template <class TValue, int ArraySize>
using EigenArray = EigenMatrix<Eigen::Matrix<TValue,ArraySize,1>>;


template <class TValue>
using EigenVector = EigenMatrix<Eigen::Matrix<TValue,Eigen::Dynamic,1>>;


} // namespace cie::linalg


#endif
