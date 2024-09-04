#ifndef CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP
#define CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ProjectiveTransformDerivative<TValue,Dimension>::evaluate(ConstIterator itBegin,
                                                          ConstIterator itEnd,
                                                          Iterator itOut) const
{
    // nD (but more importantly 3D) projective transform derivatives
    // will have to wait until I figure out how to implement them.
    static_assert(Dimension == 2, "Projective transformations are only supported in 2D for now.");

    CIE_OUT_OF_RANGE_CHECK(std::distance(itBegin, itEnd) == Dimension)
    StaticArray<TValue,Dimension+1> homogeneousInput;
    std::copy(itBegin, itEnd, homogeneousInput.begin());
    homogeneousInput.back() = 1;

    // Linear part: compute the product of the stored 3D matrix with
    // +---+---+-----+---+
    // | x   x   ...   x |
    // | y   y   ...   y |
    // | .   .   ...   . |
    // | .   .   ...   . |
    // | .   .   ...   . |
    // | 1   1   ...   1 |
    // +---+---+-----+---+
    StaticArray<TValue, Dimension*Dimension> output;
    {
        auto it = output.begin();
        typename StaticArray<TValue, Dimension>::const_iterator itInput;
        for (unsigned iComponentBegin=0; iComponentBegin<Dimension*Dimension*(Dimension+1); iComponentBegin+=(Dimension+1)) {
            itInput = homogeneousInput.begin();
            *it = 0;
            for (unsigned iDim=0; iDim<(Dimension+1); ++iDim) {
                *it += *itInput * _enumeratorCoefficients[iComponentBegin + iDim];
            } // for iDim in range(Dimension + 1)
            ++it;
        } // for iComponent in range(Dimension * Dimension * (Dimension+1), Dimension + 1)
    }

    // Nonlinear part
    const TValue denominator =   _denominatorCoefficients[0] * homogeneousInput[0]
                               + _denominatorCoefficients[1] * homogeneousInput[1]
                               + _denominatorCoefficients[2] /* * homogeneousInput[2] */;
    const TValue scale = static_cast<TValue>(1) / (denominator * denominator);

    // Scale and output
    for (unsigned i=0; i<Dimension*Dimension; ++i) {
        *itOut++ = scale * output[i];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
ProjectiveTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator itBegin,
                                                                     ConstIterator itEnd) const
{
    StaticArray<TValue,Dimension*Dimension> derivative;

    CIE_BEGIN_EXCEPTION_TRACING
    this->evaluate(itBegin, itEnd, derivative.data());
    CIE_END_EXCEPTION_TRACING

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
