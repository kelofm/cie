// --- External Includes ---
#include "Eigen/LU"

// --- FEM Includes ---
#include "packages/utilities/inc/template_macros.hpp"
#include "packages/maths/inc/AffineTransform.hpp"

// --- Linalg Includes ---
#include "packages/overloads/inc/matrix_operators.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>
#include <optional>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
AffineTransformDerivative<TValue,Dimension>::AffineTransformDerivative() noexcept
    : _matrix(TransformationMatrix::makeIdentityMatrix())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
AffineTransformDerivative<TValue,Dimension>::AffineTransformDerivative(Ref<const AffineTransform<TValue,Dimension>> rTransform)
{
    // First Dimension x Dimension sub of the original transformation matrix
    this->_matrix.wrapped() = rTransform.getTransformationMatrix().wrapped().template block<Dimension,Dimension>(0, 0);
}


namespace detail {


template <class TValue, unsigned Dimension>
class AffineCoefficients
{
public:
    using Matrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1,Dimension+1>;

public:
    AffineCoefficients()
        : _matrix()
    {
        CIE_BEGIN_EXCEPTION_TRACING

        constexpr unsigned numberOfPoints = Dimension + 1;

        // Initial point [[-1]^D, 1]
        StaticArray<TValue,Dimension+1> point;
        std::fill(point.begin(),
                  point.begin() + Dimension,
                  -1);
        point.back() = 1; // <== augmented component

        for (unsigned iPoint=0; iPoint<numberOfPoints; iPoint++) {
            // Copy point components into the coefficient matrix
            for (unsigned iComponent=0; iComponent<Dimension+1; iComponent++) {
                _matrix(iComponent, iPoint) = point[iComponent];
            }

            // Update point for the next iteration
            point[iPoint] = 1;
            if (iPoint) {
                point[iPoint - 1] = -1;
            }
        } // for iPoint in range(Dimension + 1)

        this->_matrix.wrapped() = this->_matrix.wrapped().inverse().eval();

        CIE_END_EXCEPTION_TRACING
    }

    Ref<const Matrix> get() const noexcept
    {
        return this->_matrix;
    }

private:
    Matrix _matrix;
}; // class AffineCoefficients


template <class TValue, unsigned Dimension>
class AffineCoefficientsSingleton
{
public:
    static Ref<const AffineCoefficients<TValue,Dimension>> get() noexcept
    {
        if (!_object.has_value()) {
            _object.emplace();
        }
        return _object.value();
    }

    static void clear() noexcept
    {
        _object.reset();
    }

private:
    static std::optional<AffineCoefficients<TValue,Dimension>> _object;
}; // class AffineCoefficientsSingleton


template <class TValue, unsigned Dimension>
std::optional<AffineCoefficients<TValue,Dimension>>
AffineCoefficientsSingleton<TValue,Dimension>::_object;


} // namespace detail


template <concepts::Numeric TValue, unsigned Dimension>
AffineTransform<TValue,Dimension>::AffineTransform() noexcept
    : AffineTransform(TransformationMatrix::makeIdentityMatrix())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
AffineTransform<TValue,Dimension>::AffineTransform(std::span<const Point> transformed)
    : AffineTransform()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(transformed.size() == Dimension + 1)

    // Assemble RHS
    StaticArray<TValue,(Dimension+1)*(Dimension+1)> homogeneousPoints;

    auto itTransformedBegin = transformed.begin();
    const auto itTransformedEnd = transformed.end();

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
AffineTransform<TValue,Dimension>::AffineTransform(RightRef<TransformationMatrix> rMatrix) noexcept
    : _transformationMatrix(std::move(rMatrix))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
typename AffineTransform<TValue,Dimension>::Derivative
AffineTransform<TValue,Dimension>::makeDerivative() const noexcept
{
    return AffineTransformDerivative<TValue,Dimension>(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
void
AffineTransform<TValue,Dimension>::computeTransformationMatrix(Ptr<const TValue> pTransformedBegin,
                                                               Ref<TransformationMatrix> rMatrix)
{
    CIE_BEGIN_EXCEPTION_TRACING
    Eigen::Map<const Eigen::Matrix<TValue,Dimension+1,Dimension+1>> homogeneousPoints(pTransformedBegin);
    rMatrix.wrapped().noalias() = homogeneousPoints * detail::AffineCoefficientsSingleton<TValue,Dimension>::get().get().wrapped();
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
typename AffineTransform<TValue,Dimension>::Inverse
AffineTransform<TValue,Dimension>::makeInverse() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return AffineTransform<TValue,Dimension>(
        typename AffineTransform<TValue,Dimension>::TransformationMatrix(this->getTransformationMatrix().wrapped().inverse())
    );
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
inline Ref<const typename AffineTransform<TValue,Dimension>::TransformationMatrix>
AffineTransform<TValue,Dimension>::getTransformationMatrix() const noexcept
{
    return _transformationMatrix;
}


template <concepts::Numeric TValue, unsigned Dimension>
inline Ref<typename AffineTransform<TValue,Dimension>::TransformationMatrix>
AffineTransform<TValue,Dimension>::getTransformationMatrix() noexcept
{
    return _transformationMatrix;
}


CIE_FEM_INSTANTIATE_TEMPLATE(AffineTransformDerivative);


CIE_FEM_INSTANTIATE_TEMPLATE(AffineTransform);


} // namespace cie::fem::maths
