#ifndef CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP
#define CIE_FEM_PROJECTIVE_TRANSFORM_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ProjectiveTransformDerivative<TValue,Dimension>::evaluate(ConstIterator it_begin,
                                                          ConstIterator it_end,
                                                          Iterator it_out) const
{
    // nD (but more importantly 3D) projective transform derivatives
    // will have to wait until I figure out how to implement them.
    static_assert(Dimension == 2, "Projective transformations are only supported in 2D for now.");

    CIE_OUT_OF_RANGE_CHECK(std::distance(it_begin, it_end) == Dimension)
    StaticArray<TValue,Dimension+1> homogeneousInput;
    std::copy(it_begin, it_end, homogeneousInput.begin());
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
        typename StaticArray<TValue, Dimension>::const_iterator it_input;
        for (unsigned i_componentBegin=0; i_componentBegin<Dimension*Dimension*(Dimension+1); i_componentBegin+=(Dimension+1)) {
            it_input = homogeneousInput.begin();
            *it = 0;
            for (unsigned i_dim=0; i_dim<(Dimension+1); ++i_dim) {
                *it += *it_input * _enumeratorCoefficients[i_componentBegin + i_dim];
            } // for i_dim in range(Dimension + 1)
            ++it;
        } // for i_component in range(Dimension * Dimension * (Dimension+1), Dimension + 1)
    }

    // Nonlinear part
    const TValue denominator =   _denominatorCoefficients[0] * homogeneousInput[0]
                               + _denominatorCoefficients[1] * homogeneousInput[1]
                               + _denominatorCoefficients[2] /* * homogeneousInput[2] */;
    const TValue scale = static_cast<TValue>(1) / (denominator * denominator);

    // Scale and output
    for (unsigned i=0; i<Dimension*Dimension; ++i) {
        *it_out++ = scale * output[i];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
inline TValue
ProjectiveTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstIterator it_begin,
                                                                     ConstIterator it_end) const
{
    StaticArray<TValue,Dimension*Dimension> derivative;

    CIE_BEGIN_EXCEPTION_TRACING
    this->evaluate(it_begin, it_end, derivative.data());
    CIE_END_EXCEPTION_TRACING

    return Eigen::Map<Eigen::Matrix<TValue,Dimension,Dimension>>(derivative.data()).determinant();
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform(TPointIt it_transformedBegin,
                                                           TPointIt it_transformedEnd)
    : ProjectiveTransform()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(std::distance(it_transformedBegin, it_transformedEnd) == Dimension*Dimension)

    // Assemble RHS
    StaticArray<TValue,Dimension*Dimension*(Dimension+1)> homogeneousPoints;

    // Copy transformed components to the first {{Dimension}} rows
    for (Size i_point=0 ; it_transformedBegin!=it_transformedEnd; it_transformedBegin++, i_point++) {
        CIE_OUT_OF_RANGE_CHECK(Dimension <= it_transformedBegin->size())
        const auto i_componentBegin = i_point * (Dimension + 1);
        for (Size i_component=0; i_component<Dimension; i_component++) {
            // This array will be interpreted as an eigen matrix, which
            // stores its data columnwise by default, so the order of the
            // components must follow that.
            homogeneousPoints[i_componentBegin + i_component] = it_transformedBegin->at(i_component);
        } // for component in point
        homogeneousPoints[i_componentBegin + Dimension] = 1; // <== last row contains homogeneous components
    } // for point in transformedPoints

    // Solve for transformation matrix components
    this->computeTransformationMatrix(homogeneousPoints.data(),
                                      this->getTransformationMatrix());

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
inline void
ProjectiveTransform<TValue,Dimension>::evaluate(ConstIterator it_argumentBegin,
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

    // Dehomogenize
    const TValue homogeneousComponent = transformed[Dimension];
    CIE_DIVISION_BY_ZERO_CHECK(homogeneousComponent != 0)
    const TValue scale = static_cast<TValue>(1) / homogeneousComponent;

    std::transform(
        transformed.begin(),
        transformed.begin() + Dimension,
        it_out,
        [scale](TValue component) {return component * scale;}
    );
}


} // namespace cie::fem::maths


#endif
