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
ProjectiveTransformDerivative<TValue,Dimension>::ProjectiveTransformDerivative() noexcept
    : _enumeratorCoefficients(),
      _denominatorCoefficients({0, 0, 1})
{
    // Set w components on the main diagonal to unity, everything else to zero
    std::fill(_enumeratorCoefficients.begin(),
              _enumeratorCoefficients.end(),
              0);

    for (unsigned i_component=Dimension; i_component<_enumeratorCoefficients.size(); i_component+=(Dimension+1)*Dimension) {
        _enumeratorCoefficients[i_component] = 1;
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransformDerivative<TValue,Dimension>::ProjectiveTransformDerivative(Ref<const ProjectiveTransform<TValue,Dimension>> r_projection)
    : _enumeratorCoefficients(),
      _denominatorCoefficients()
{
    static_assert(Dimension == 2, "Projective transformations are only supported in 2D for now.");

    // Compute denominator coefficients

    const auto& r_matrix = r_projection.getTransformationMatrix();
    for (unsigned i_coefficient=0; i_coefficient<=Dimension; ++i_coefficient) {
        _denominatorCoefficients[i_coefficient] = r_matrix(Dimension, i_coefficient);
    }

    // Suppose the the projective transform matrix looks like this:
    // (ignore the discrepancy between the row-major naming here
    // and the general column-wise data storage in the implementation).
    // +---+---+---+
    // | a   b   c |
    // | d   e   f |
    // | g   h   i |
    // +---+---+---+
    // The derivative's 3D enumerator matrix then follows from subdeterminants:
    // + ----------------------+-----------------------+
    // | [    0, ah-bg,   i-c]   [    0, dh-eg,   i-f] |
    // | [bg-ah,     0,   i-c]   [eg-dh,     0,   i-f] |
    // +-----------------------+-----------------------+

    // Compute temporaries

    // ah-bg
    const TValue ahbg =   r_matrix(0, 0) * r_matrix(Dimension, 1)
                        - r_matrix(0, 1) * r_matrix(Dimension, 0);

    // dh-eg
    const TValue dheg =   r_matrix(1, 0) * r_matrix(Dimension, 1)
                        - r_matrix(1, 1) * r_matrix(Dimension, 0);

    // i-c
    const TValue omc = r_matrix(Dimension, Dimension) - r_matrix(0, Dimension);

    // i-f
    const TValue omf = r_matrix(Dimension, Dimension) - r_matrix(1, Dimension);

    // Compute enumerator coefficients

    // [0, 0, :]
    _enumeratorCoefficients[ 0] = 0;
    _enumeratorCoefficients[ 1] = ahbg;
    _enumeratorCoefficients[ 2] = omc;

    // [1, 0, :]
    _enumeratorCoefficients[ 3] = -ahbg;
    _enumeratorCoefficients[ 4] = 0;
    _enumeratorCoefficients[ 5] = omc;

    // [0, 1, :]
    _enumeratorCoefficients[ 6] = 0;
    _enumeratorCoefficients[ 7] = dheg;
    _enumeratorCoefficients[ 8] = omf;

    // [1, 1, :]
    _enumeratorCoefficients[ 9] = -dheg;
    _enumeratorCoefficients[10] = 0;
    _enumeratorCoefficients[11] = omf;
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

        for (unsigned i_point=0; i_point<Dimension+1; i_point++, ++permutation) {
            for (unsigned i_component=0; i_component<Dimension; i_component++) {
                _matrix(i_component, i_point) = *(*permutation)[i_component];
            }
            _matrix(Dimension, i_point) = 1;
        }

        // Right hand side == [1]^(D+1)
        Eigen::Matrix<TValue,Dimension+1,1> rhs;
        for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
            rhs(i_dim, 0) = *(*permutation)[i_dim];
        }
        rhs(Dimension, 0) = 1;

        // Solve for column coefficients and scale the columns
        Eigen::Matrix<TValue,Dimension+1,1> columnCoefficients = _matrix.wrapped().fullPivLu().solve(rhs);
        for (unsigned i_row=0; i_row<Dimension+1; ++i_row) {
            for (unsigned i_column=0; i_column<Dimension+1; ++i_column) {
                _matrix(i_row, i_column) *= columnCoefficients(i_column, 0);
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


template <class TValue, unsigned Dimension>
struct ComputeProjectiveMatrix
{
    static void compute(Ptr<TValue> p_transformedBegin,
                        Ref<typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix> r_matrix)
    {
        throw NotImplementedException("","");
    }
};


template <class TValue>
struct ComputeProjectiveMatrix<TValue,2>
{
    static void compute(Ptr<TValue> p_transformedBegin,
                        Ref<typename ProjectiveTransform<TValue,2>::TransformationMatrix> r_matrix)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        constexpr unsigned Dimension = 2;
        Eigen::Map<Eigen::Matrix<TValue,Dimension+1,Dimension+1>> homogeneousPoints(p_transformedBegin);
        Eigen::Map<const Eigen::Matrix<TValue,Dimension+1,1>> rhs(p_transformedBegin + (Dimension + 1) * (Dimension + 1));

        const Eigen::Matrix<TValue,Dimension+1,1> homogeneousSolution = Eigen::FullPivLU<Eigen::Matrix<TValue,Dimension+1,Dimension+1>>(homogeneousPoints).solve(rhs);
        for (unsigned i_point=0; i_point<Dimension+1; i_point++) {
            const TValue scale = homogeneousSolution[i_point];
            for (unsigned i_component=0; i_component<Dimension+1; i_component++) {
                homogeneousPoints(i_component, i_point) *= scale;
            }
        }

        r_matrix.wrapped().noalias() = homogeneousPoints * detail::ProjectiveCoefficientsSingleton<TValue,Dimension>::get().get().wrapped();
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
ProjectiveTransform<TValue,Dimension>::ProjectiveTransform(RightRef<TransformationMatrix> r_matrix) noexcept
    : _transformationMatrix(std::move(r_matrix))
{
}


template <concepts::Numeric TValue, unsigned Dimension>
void
ProjectiveTransform<TValue,Dimension>::computeTransformationMatrix(Ptr<TValue> p_transformedBegin,
                                                                   Ref<TransformationMatrix> r_matrix)
{
    CIE_BEGIN_EXCEPTION_TRACING
    detail::ComputeProjectiveMatrix<TValue,Dimension>::compute(p_transformedBegin, r_matrix);
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ProjectiveTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransform<TValue,Dimension>
ProjectiveTransform<TValue,Dimension>::makeInverse() const
{
    CIE_BEGIN_EXCEPTION_TRACING
    return ProjectiveTransform<TValue,Dimension>(
        typename ProjectiveTransform<TValue,Dimension>::TransformationMatrix(this->getTransformationMatrix().wrapped().inverse())
    );
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
ProjectiveTransformDerivative<TValue,Dimension>
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


template class ProjectiveTransformDerivative<float,2>;


template class ProjectiveTransformDerivative<double,2>;


template class ProjectiveTransform<float,2>;


template class ProjectiveTransform<double,2>;


} // namespace cie::fem::maths
