#ifndef CIE_LINALG_UTILITIES_HPP
#define CIE_LINALG_UTILITIES_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- Internal includes ---
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/utilities/inc/errortypes.hpp"

namespace cie::linalg {


using PermutationVector = std::vector<Size>;

void updatePermutation(const DynamicEigenMatrix<double>& matrix,
                       PermutationVector& permutation,
                       Size index,
                       double singularTolerance);


} // namespace cie::linalg


#endif
