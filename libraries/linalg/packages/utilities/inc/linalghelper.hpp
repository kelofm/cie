#ifndef CIE_LINALG_HELPER_HPP
#define CIE_LINALG_HELPER_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- Internal Includes ---
#include "packages/types/inc/vectortypes.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"

// --- STL Includes ---
#include <iostream>
#include <stdexcept>


namespace cie::linalg::linalghelper {


template<class RowFunction>
void writeRow(const RowFunction& vector, Size size, std::ostream& out, Size digits);

template<concepts::NumericContainer ContainerType>
void write(const ContainerType& vector, std::ostream& out = std::cout);

void write(const DynamicEigenMatrix<double>& matrix, std::ostream& out = std::cout);


} // namespace cie::linalg::linalghelper

#include "packages/utilities/impl/linalghelper_impl.hpp"

#endif
