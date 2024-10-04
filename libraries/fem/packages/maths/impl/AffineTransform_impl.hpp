#ifndef CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_IMPL_HPP
#define CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_IMPL_HPP

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
AffineTransformDerivative<TValue,Dimension>::evaluate(ConstIterator,
                                                      ConstIterator,
                                                      Iterator itOut) const
{
    std::copy(this->_matrix.wrapped().data(),
              this->_matrix.wrapped().data() + Dimension * Dimension,
              itOut);
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
AffineTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator,
                                                                 ConstIterator) const
{
    return this->_matrix.wrapped().determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator PointIterator>
AffineTransform<TValue,Dimension>::AffineTransform(PointIterator itTransformedBegin,
                                                   PointIterator itTransformedEnd)
    : AffineTransform()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(std::distance(itTransformedBegin, itTransformedEnd) == Dimension + 1)

    // Assemble RHS
    StaticArray<TValue,(Dimension+1)*(Dimension+1)> homogeneousPoints;

    // Copy transformed components to the first {{Dimension}} rows
    for (Size iPoint=0 ; itTransformedBegin!=itTransformedEnd; itTransformedBegin++, iPoint++) {
        CIE_OUT_OF_RANGE_CHECK(Dimension <= itTransformedBegin->size())
        for (Size iComponent=0; iComponent<Dimension; iComponent++) {
            // This array will be interpreted as an eigen matrix, which
            // stores its data columnwise by default, so the order of the
            // components must follow that.
            homogeneousPoints[iComponent + iPoint * (Dimension + 1)] = itTransformedBegin->at(iComponent);
        } // for component in point
        homogeneousPoints[Dimension + iPoint * (Dimension + 1)] = 1; // <== last row contains homogeneous components
    } // for point in transformedPoints

    // Solve for transformation matrix components
    this->computeTransformationMatrix(homogeneousPoints.data(),
                                      this->getTransformationMatrix());

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
AffineTransform<TValue,Dimension>::evaluate(ConstIterator itArgumentBegin,
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

    // Output result components
    std::copy(transformed.begin(),
              transformed.begin() + Dimension,
              itOut);
}


} // namespace cie::fem::maths


#endif
