#ifndef CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_IMPL_HPP
#define CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/AffineTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
AffineTransformDerivative<TValue,Dimension>::evaluate(ConstIterator it_begin,
                                                      ConstIterator it_end,
                                                      Iterator it_out) const
{
    std::copy(this->_matrix.wrapped().data(),
              this->_matrix.wrapped().data() + Dimension * Dimension,
              it_out);
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
AffineTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator it_begin,
                                                                 ConstIterator it_end) const
{
    return this->_matrix.wrapped().determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator PointIterator>
AffineTransform<TValue,Dimension>::AffineTransform(PointIterator it_transformedBegin,
                                                   PointIterator it_transformedEnd)
    : AffineTransform()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == Dimension + 1)

    // Assemble RHS
    StaticArray<TValue,(Dimension+1)*(Dimension+1)> homogeneousPoints;

    // Copy transformed components to the first {{Dimension}} rows
    for (Size i_point=0 ; it_transformedBegin!=it_transformedEnd; it_transformedBegin++, i_point++) {
        CIE_OUT_OF_RANGE_CHECK(Dimension <= it_transformedBegin->size())
        for (Size i_component=0; i_component<Dimension; i_component++) {
            // This array will be interpreted as an eigen matrix, which
            // stores its data columnwise by default, so the order of the
            // components must follow that.
            homogeneousPoints[i_component + i_point * (Dimension + 1)] = it_transformedBegin->at(i_component);
        } // for component in point
        homogeneousPoints[Dimension + i_point * (Dimension + 1)] = 1; // <== last row contains homogeneous components
    } // for point in transformedPoints

    // Solve for transformation matrix components
    this->computeTransformationMatrix(homogeneousPoints.data(),
                                      this->getTransformationMatrix());

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
AffineTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                            ConstIterator it_argumentEnd,
                                            Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(Dimension == std::distance(it_argumentBegin, it_argumentEnd))

    // Copy augmented point
    typename Kernel<Dimension,TValue>::template static_array<Dimension+1> augmentedPoint;
    std::copy(it_argumentBegin,
              it_argumentEnd,
              augmentedPoint.begin());

    augmentedPoint[Dimension] = static_cast<TValue>(1);

    // Transform
    const auto transformed = this->getTransformationMatrix() * augmentedPoint;

    // Output result components
    std::copy(
        transformed.begin(),
        transformed.begin() + Dimension,
        it_out
    );
}


} // namespace cie::fem::maths


#endif
