#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- STL Includes ---
#include <numeric> // std::inner_product


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
void ProjectiveTransformDerivative<TValue,Dimension>::evaluate(ConstSpan in, Span out) const
{
    // Transform the input point to homogenized space.
    StaticArray<TValue,Dimension+1> homogenizedInput;
    std::copy(in.begin(), in.end(), homogenizedInput.begin());
    homogenizedInput.back() = static_cast<TValue>(1);

    // Compute the denominator shared by all entries in the output matrix.
    const auto lastRow = _projectionMatrix.wrapped()(Dimension, Eigen::placeholders::all);
    const TValue lastRowProduct = std::inner_product(homogenizedInput.begin(),
                                                     homogenizedInput.end(),
                                                     lastRow.begin(),
                                                     static_cast<TValue>(0));
    const TValue denominator = std::pow(lastRowProduct, static_cast<TValue>(2));

    CIE_DIVISION_BY_ZERO_CHECK(denominator)
    const TValue scale = static_cast<TValue>(1) / denominator;

    auto itOut = out.data();
    for (unsigned iRow=0u; iRow<Dimension; ++iRow) {
        const auto row = _projectionMatrix.wrapped()(iRow, Eigen::placeholders::all);
        const TValue rowProduct = std::inner_product(
            homogenizedInput.begin(),
            homogenizedInput.end(),
            row.begin(),
            static_cast<TValue>(0));

        for (unsigned iColumn=0u; iColumn<Dimension; ++iColumn) {
            *itOut++ = scale * (_projectionMatrix.wrapped()(iRow, iColumn) * lastRowProduct - _projectionMatrix(Dimension, iColumn) * rowProduct);
        } // for iColumn in range(Dimension)
    } // for iRow in range(Dimension)
}


template <concepts::Numeric TValue, unsigned Dimension>
TValue ProjectiveTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstSpan in) const
{
    StaticArray<TValue,Dimension*Dimension> derivative;
    this->evaluate(in, derivative);
    return Eigen::Map<Eigen::Matrix<TValue,Dimension,Dimension>>(derivative.data()).determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
constexpr unsigned ProjectiveTransformDerivative<TValue,Dimension>::size() noexcept {
    return Dimension * Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
constexpr unsigned ProjectiveTransform<TValue,Dimension>::size() noexcept {
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
void ProjectiveTransform<TValue,Dimension>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(Dimension == in.size())
    CIE_OUT_OF_RANGE_CHECK(Dimension == out.size())

    // Copy augmented point
    typename Kernel<Dimension,TValue>::template static_array<Dimension+1> augmentedPoint;
    for (Size iDim=0; iDim<Dimension; ++iDim) {
        augmentedPoint[iDim] = in[iDim];
    }

    // <== GCC thinks this doesn't initialize augmentedPoint ...
    //std::copy(itArgumentBegin,
    //          itArgumentEnd,
    //          augmentedPoint.begin());

    augmentedPoint[Dimension] = static_cast<TValue>(1);

    // Transform
    const auto transformed = this->getTransformationMatrix() * augmentedPoint;

    // Dehomogenize
    const TValue homogeneousComponent = transformed[Dimension];
    CIE_DIVISION_BY_ZERO_CHECK(homogeneousComponent != 0)
    const TValue scale = static_cast<TValue>(1) / homogeneousComponent;

    std::transform(
        transformed.begin(),
        transformed.begin() + Dimension,
        out.data(),
        [scale](TValue component) {return component * scale;}
    );
}


} // namespace cie::fem::maths
