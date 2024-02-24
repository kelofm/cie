#ifndef CIE_FEM_PROJECTIVE_TRANSFORM_HPP
#define CIE_FEM_PROJECTIVE_TRANSFORM_HPP

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransform;


/// @addtogroup fem
/// @{


/** @brief Expression representing the derivative of @ref ProjectiveTransform.
 *  @todo Implement in 3D (or any dimension).
 *  @details The derivative is a @f$ D x D @f$ matrix that can be separated into
 *           an enumerator matrix, which is a matrix product of a 3D coefficient
 *           matrix and a 2D square matrix consisting of repeated copies of the
 *           homogeneous coordinates, and a scalar denominator, which is equal
 *           to the square of the inner product of the projective transform's
 *           last row and the homogeneous coordinates.
 *           @f[
 *           \frac{\mathbf C \begin{bmatrix} \mathbf {\hat{x}} & \mathbf {\hat{x}} & \ldots & \mathbf {\hat{x}} \end{bmatrix}}
 *                {\begin{bmatrix} last\_row\_of\_ \mathbf T \end{bmatrix}
 *                 \begin{bmatrix}
 *                      x \\
 *                      y \\
 *                 \vdots \\
 *                      1
 *                 \end{bmatrix}}
 *           @f]
 *  @details The 3D coefficient matrix is stored as a flattened matrix in the
 *           following index order
 *           1) homogeneous component index (x, y, z, ..., w)
 *           2) row index (Eigen stores column-major matrices by default)
 *           3) column index
 *  @details The 2D matrix of homogeneous coordinates takes the following form
 *           @f[ \begin{bmatrix} \mathbf {\hat{x}} & \mathbf {\hat{x}} & \ldots & \mathbf {\hat{x}} \end{bmatrix} = \begin{bmatrix}
 *                x &      x & \ldots &      x \\
 *                y &      y & \ldots &      y \\
 *           \vdots & \vdots & \ddots & \vdots \\
 *                1 &      1 & \ldots &      1
 *           \end{bmatrix} @f]
 *
 */
template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransformDerivative : public ExpressionTraits<TValue>
{
private:
    /** @brief Flattened 3D matrix storing the enumerator's coefficients.
     *  @details The enumerator of the derivative is computed by evaluating the product
     *           of this 3D matrix with the following 2D matrix:
     *           @code
     *           +---+---+-----+---+
     *           | x   x   ...   x |
     *           | y   y   ...   y |
     *           | .   .   ...   . |
     *           | .   .   ...   . |
     *           | .   .   ...   . |
     *           | 1   1   ...   1 |
     *           +---+---+-----+---+
     *           @endcode
     *  @details Finally, we can get the derivative by dividing each component of the
     *           product computed in the previous step by the following denominator:
     *           @code
     *                                          +---+
     *           [_denominatorCoefficients]  *  | x |
     *                                          | y |
     *                                          | . |
     *                                          | . |
     *                                          | . |
     *                                          | 1 |
     *                                          +---+
     *           @endcode
     */
    using EnumeratorCoefficients = StaticArray<TValue,Dimension*Dimension*(Dimension+1)>;

public:
    CIE_DEFINE_CLASS_POINTERS(ProjectiveTransformDerivative)

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity by default.
    ProjectiveTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(ConstIterator it_begin,
                  ConstIterator it_end,
                  Iterator it_out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Compute the determinant of the projective transform's jacobian.
    TValue evaluateDeterminant(ConstIterator it_begin, ConstIterator it_end) const;

private:
    friend class ProjectiveTransform<TValue,Dimension>;

    /// @brief Construct from a @ref ProjectiveTransform.
    ProjectiveTransformDerivative(Ref<const ProjectiveTransform<TValue,Dimension>> r_projection);

private:
    EnumeratorCoefficients _enumeratorCoefficients;

    StaticArray<TValue,3> _denominatorCoefficients;
}; // class ProjectiveTransformDerivative



/** @brief Class representing a projective transform.
 *  @details Uniquely defines a mapping between any pair of generalized quads.
 *  @note Unfortunately, I haven't found the proper term for quadrilaterals
 *        generalized to arbitrary dimensions, so @a generalized @a quad
 *        will have to do for now.
 */
template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransform : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(ProjectiveTransform)

    using TransformationMatrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1, Dimension+1>;

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity transform by default
    ProjectiveTransform() noexcept;

    /** @brief Projective transformation from @f$ 2 \cdot D @f$ transformed points.
     *
     *  @details The transformation is characterized by how it deforms a cube
     *           defined on @f$ [-1,1]^D @f$.
     *
     *  @details To avoid unintended overlapping transformations, the order of vertices
     *           should match an outer product with smaller spatial indices varying more
     *           frequently. Example in 3D:
     *           @f[ \begin{bmatrix}
     *           -1 &  1 & -1 &  1 & -1 &  1 & -1 &  1 \\
     *           -1 & -1 &  1 &  1 & -1 & -1 &  1 &  1 \\
     *           -1 & -1 & -1 & -1 &  1 &  1 &  1 &  1
     *           \end{bmatrix} @f]
     *
     *  @param it_transformedBegin iterator pointing to the transformed cube's base @f$ [-1]^D @f$.
     *  @param it_transformedEnd iterator past the last transformed point (should be identical to it_transformedBegin + 2 * D + 1).
     */
    template <concepts::Iterator PointIt>
    ProjectiveTransform(PointIt it_transformedBegin,
                        PointIt it_transformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components.
    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the inverse transform.
    ProjectiveTransform makeInverse() const;

    /// @brief Construct the derivative of the projective transform.
    ProjectiveTransformDerivative<TValue,Dimension> makeDerivative() const;

    /// @brief Get the matrix representation of the transformation.
    Ref<const TransformationMatrix> getTransformationMatrix() const noexcept;

private:
    /// @brief Construct from a precomputed transformation matrix.
    ProjectiveTransform(RightRef<TransformationMatrix> r_matrix) noexcept;

    /// @brief Get the matrix representation of the transformation.
    Ref<TransformationMatrix> getTransformationMatrix() noexcept;

    /// @brief Compute the transformation matrix from the homogeneous representation of transformed points.
    /// @param[in] p_transformedBegin Ptr to the first component of the homogenized transformed points.
    /// @param[out] r_matrix Transformation matrix to write to.
    /// @warning @a p_transformedBegin is mutated during construction.
    static void computeTransformationMatrix(Ptr<TValue> p_transformedBegin,
                                            Ref<TransformationMatrix> r_matrix);

private:
    TransformationMatrix _transformationMatrix;
}; // class ProjectiveTransform


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/ProjectiveTransform_impl.hpp"

#endif
