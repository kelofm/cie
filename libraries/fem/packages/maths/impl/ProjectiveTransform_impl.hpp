#ifndef CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP
#define CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- STL Includes ---
#include <numeric> // std::inner_product


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ProjectiveTransformDerivative<TValue,Dimension>::evaluate(ConstIterator itBegin,
                                                          ConstIterator itEnd,
                                                          Iterator itOut) const
{
    // Transform the input point to homogenized space.
    StaticArray<TValue,Dimension+1> homogenizedInput;
    std::copy(itBegin, itEnd, homogenizedInput.begin());
    homogenizedInput.back() = static_cast<TValue>(1);

    // Compute the denominator shared by all entries in the output matrix.
    const auto lastRow = _projectionMatrix.wrapped()(Dimension, Eigen::all);
    const TValue lastRowProduct = std::inner_product(homogenizedInput.begin(),
                                                     homogenizedInput.end(),
                                                     lastRow.begin(),
                                                     static_cast<TValue>(0));
    const TValue denominator = std::pow(lastRowProduct, static_cast<TValue>(2));

    CIE_DIVISION_BY_ZERO_CHECK(denominator)
    const TValue scale = static_cast<TValue>(1) / denominator;

    for (unsigned iRow=0u; iRow<Dimension; ++iRow) {
        const auto row = _projectionMatrix.wrapped()(iRow, Eigen::all);
        const TValue rowProduct = std::inner_product(
            homogenizedInput.begin(),
            homogenizedInput.end(),
            row.begin(),
            static_cast<TValue>(0));

        for (unsigned iColumn=0u; iColumn<Dimension; ++iColumn) {
            *itOut++ = scale * (_projectionMatrix.wrapped()(iRow, iColumn) * lastRowProduct - rowProduct);
        } // for iColumn in range(Dimension)
    } // for iRow in range(Dimension)
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
ProjectiveTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator itBegin,
                                                                     ConstIterator itEnd) const
{
    StaticArray<TValue,Dimension*Dimension> derivative;
    this->evaluate(itBegin, itEnd, derivative.data());
    return Eigen::Map<Eigen::Matrix<TValue,Dimension,Dimension>>(derivative.data()).determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform(TPointIt itTransformedBegin,
                                                           TPointIt itTransformedEnd)
    : ProjectiveTransform()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(std::distance(itTransformedBegin, itTransformedEnd) == Dimension*Dimension)

    // Assemble RHS
    StaticArray<TValue,Dimension*Dimension*(Dimension+1)> homogeneousPoints;

    // Copy transformed components to the first {{Dimension}} rows
    for (Size iPoint=0 ; itTransformedBegin!=itTransformedEnd; itTransformedBegin++, iPoint++) {
        CIE_OUT_OF_RANGE_CHECK(Dimension <= itTransformedBegin->size())
        const auto iComponentBegin = iPoint * (Dimension + 1);
        for (Size iComponent=0; iComponent<Dimension; iComponent++) {
            // This array will be interpreted as an eigen matrix, which
            // stores its data columnwise by default, so the order of the
            // components must follow that.
            homogeneousPoints[iComponentBegin + iComponent] = itTransformedBegin->at(iComponent);
        } // for component in point
        homogeneousPoints[iComponentBegin + Dimension] = 1; // <== last row contains homogeneous components
    } // for point in transformedPoints

    // Solve for transformation matrix components
    this->computeTransformationMatrix(homogeneousPoints.data(),
                                      this->getTransformationMatrix());

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ProjectiveTransform<TValue,Dimension>::evaluate(ConstIterator itArgumentBegin,
                                                [[maybe_unused]] ConstIterator itArgumentEnd,
                                                Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(Dimension == std::distance(itArgumentBegin, itArgumentEnd))

    // Copy augmented point
    typename Kernel<Dimension,TValue>::template static_array<Dimension+1> augmentedPoint;
    for (Size iDim=0; iDim<Dimension; ++iDim) {
        augmentedPoint[iDim] = itArgumentBegin[iDim];
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
        itOut,
        [scale](TValue component) {return component * scale;}
    );
}


} // namespace cie::fem::maths


#endif
