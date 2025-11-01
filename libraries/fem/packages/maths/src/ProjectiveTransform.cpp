// --- External Includes ---
#include "Eigen/LU"

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/stl_extension/inc/state_iterator.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <optional>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform(std::span<const Point> transformed)
    : ProjectiveTransform()
{
    CIE_OUT_OF_RANGE_CHECK(transformed.size() == Dimension * Dimension)

    CIE_BEGIN_EXCEPTION_TRACING

    // Assemble RHS
    auto itTransformedBegin = transformed.begin();
    const auto itTransformedEnd = transformed.end();
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
ProjectiveTransformDerivative<TValue,Dimension>::ProjectiveTransformDerivative() noexcept
    : ProjectiveTransformDerivative(TransformationMatrix::makeIdentityMatrix())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransformDerivative<TValue,Dimension>::ProjectiveTransformDerivative(Ref<const ProjectiveTransform<TValue,Dimension>> rProjection)
    : ProjectiveTransformDerivative(rProjection.getTransformationMatrix())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransformDerivative<TValue,Dimension>::ProjectiveTransformDerivative(Ref<const TransformationMatrix> rProjectionMatrix) noexcept
    : _projectionMatrix(rProjectionMatrix)
{
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ProjectiveTransformDerivative<TValue,Dimension>::size() const noexcept
{
    return Dimension * Dimension;
}


namespace detail {


template <class TValue, unsigned Dimension>
class ProjectiveCoefficients
{
public:
    using Matrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1,Dimension+1>;

public:
    ProjectiveCoefficients()
        : _matrix()
    {
        CIE_BEGIN_EXCEPTION_TRACING

        /// @todo Implement for higher dimensions (currently 2D only).
        StaticArray<TValue,2> states {-1, 1};
        auto permutation = utils::makeInternalStateIterator(states, Dimension);

        for (unsigned iPoint=0; iPoint<Dimension+1; iPoint++, ++permutation) {
            for (unsigned iComponent=0; iComponent<Dimension; iComponent++) {
                _matrix(iComponent, iPoint) = *(*permutation)[iComponent];
            }
            _matrix(Dimension, iPoint) = 1;
        }

        // Right hand side == [1]^(D+1)
        Eigen::Matrix<TValue,Dimension+1,1> rhs;
        for (unsigned iDim=0; iDim<Dimension; ++iDim) {
            rhs(iDim, 0) = *(*permutation)[iDim];
        }
        rhs(Dimension, 0) = 1;

        // Solve for column coefficients and scale the columns
        Eigen::Matrix<TValue,Dimension+1,1> columnCoefficients = _matrix.wrapped().fullPivLu().solve(rhs);
        for (unsigned iRow=0; iRow<Dimension+1; ++iRow) {
            for (unsigned iColumn=0; iColumn<Dimension+1; ++iColumn) {
                _matrix(iRow, iColumn) *= columnCoefficients(iColumn, 0);
            }
        }

        // Invert the result
        this->_matrix.wrapped() = this->_matrix.wrapped().inverse().eval();

        CIE_END_EXCEPTION_TRACING
    }

    Ref<const Matrix> get() const noexcept
    {
        return this->_matrix;
    }

private:
    Matrix _matrix;
}; // class ProjectiveCoefficients


template <class TValue, unsigned Dimension>
class ProjectiveCoefficientsSingleton
{
public:
    static Ref<const ProjectiveCoefficients<TValue,Dimension>> get() noexcept
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
    static std::optional<ProjectiveCoefficients<TValue,Dimension>> _object;
}; // class ProjectiveCoefficientsSingleton


template <class TValue, unsigned Dimension>
std::optional<ProjectiveCoefficients<TValue,Dimension>>
ProjectiveCoefficientsSingleton<TValue,Dimension>::_object;


//template <class TValue, unsigned Dimension>
//struct ComputeProjectiveMatrix
//{
//    static void compute(Ptr<TValue>, Ref<typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix>)
//    {
//        throw NotImplementedException("","");
//    }
//};


//template <class TValue>
//struct ComputeProjectiveMatrix<TValue,2>
template <class TValue, unsigned Dimension>
struct ComputeProjectiveMatrix
{
    static void compute(Ptr<TValue> pTransformedBegin,
                        Ref<typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix> rMatrix)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        //constexpr unsigned Dimension = 2;
        Eigen::Map<Eigen::Matrix<TValue,Dimension+1,Dimension+1>> homogeneousPoints(pTransformedBegin);
        Eigen::Map<const Eigen::Matrix<TValue,Dimension+1,1>> rhs(pTransformedBegin + (Dimension + 1) * (Dimension + 1));

        Eigen::FullPivLU<Eigen::Matrix<TValue,Dimension+1,Dimension+1>> solver;
        solver.compute(homogeneousPoints);

        if (!solver.isInvertible()) {
            CIE_THROW(
                Exception,
                "Singular input for projective transform.\n"
                << "LHS:\n" << homogeneousPoints << "\n"
                << "RHS:\n" << rhs)
        }

        const Eigen::Matrix<TValue,Dimension+1,1> homogeneousSolution = solver.solve(rhs);

        for (unsigned iPoint=0; iPoint<Dimension+1; iPoint++) {
            const TValue scale = homogeneousSolution[iPoint];
            for (unsigned iComponent=0; iComponent<Dimension+1; iComponent++) {
                homogeneousPoints(iComponent, iPoint) *= scale;
            }
        }

        rMatrix.wrapped().noalias() = homogeneousPoints * detail::ProjectiveCoefficientsSingleton<TValue,Dimension>::get().get().wrapped();
        CIE_END_EXCEPTION_TRACING
    }
};


} // namespace detail


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform() noexcept
    : ProjectiveTransform(TransformationMatrix::makeIdentityMatrix())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform(RightRef<TransformationMatrix> rMatrix) noexcept
    : _transformationMatrix(std::move(rMatrix))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
void
ProjectiveTransform<TValue,Dimension>::computeTransformationMatrix(Ptr<TValue> pTransformedBegin,
                                                                   Ref<TransformationMatrix> rMatrix)
{
    CIE_BEGIN_EXCEPTION_TRACING
    detail::ComputeProjectiveMatrix<TValue,Dimension>::compute(pTransformedBegin, rMatrix);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ProjectiveTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
typename ProjectiveTransform<TValue,Dimension>::Inverse
ProjectiveTransform<TValue,Dimension>::makeInverse() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return ProjectiveTransform<TValue,Dimension>(
        typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix(this->getTransformationMatrix().wrapped().inverse())
    );
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
typename ProjectiveTransform<TValue,Dimension>::Derivative
ProjectiveTransform<TValue,Dimension>::makeDerivative() const
{
    return ProjectiveTransformDerivative<TValue,Dimension>(*this);
}


template <concepts::Numeric TValue, unsigned Dimension>
inline Ref<const typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix>
ProjectiveTransform<TValue,Dimension>::getTransformationMatrix() const noexcept
{
    return _transformationMatrix;
}


template <concepts::Numeric TValue, unsigned Dimension>
inline Ref<typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix>
ProjectiveTransform<TValue,Dimension>::getTransformationMatrix() noexcept
{
    return _transformationMatrix;
}


CIE_FEM_INSTANTIATE_TEMPLATE(ProjectiveTransformDerivative);


CIE_FEM_INSTANTIATE_TEMPLATE(ProjectiveTransform);


} // namespace cie::fem::maths
