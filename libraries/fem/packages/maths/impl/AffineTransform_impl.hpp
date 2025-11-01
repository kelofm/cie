#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/AffineTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
AffineTransformDerivative<TValue,Dimension>::evaluate(ConstSpan, Span output) const
{
    std::copy(this->_matrix.wrapped().data(),
              this->_matrix.wrapped().data() + Dimension * Dimension,
              output.data());
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
AffineTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstSpan) const
{
    return this->_matrix.wrapped().determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
AffineTransform<TValue,Dimension>::evaluate(ConstSpan input, Span output) const
{
    CIE_OUT_OF_RANGE_CHECK(Dimension == input.size())
    CIE_OUT_OF_RANGE_CHECK(Dimension == output.size())

    // Copy augmented point
    typename Kernel<Dimension,TValue>::template static_array<Dimension+1> augmentedPoint;
    std::copy(input.data(),
              input.data() + Dimension,
              augmentedPoint.data());

    augmentedPoint[Dimension] = static_cast<TValue>(1);

    // Transform
    const auto transformed = this->getTransformationMatrix() * augmentedPoint;

    // Output result components
    std::copy(transformed.begin(),
              transformed.begin() + Dimension,
              output.data());
}


} // namespace cie::fem::maths
