#pragma once

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransform;


/// @addtogroup fem
/// @{


/** @brief Expression representing the derivative of @ref ProjectiveTransform.
 *  @details Given the transformation matrix of a @ref ProjectiveTransform "projective transform" @f$t_{ij}@f$,
 *           in @f$d@f$ dimensions and the input vector @f$x@f$ in homogenized space
 *           @f$\xi = \begin{bmatrix} x^T & 1 \end{bmatrix}^T,
 *           @f$ this class is responsible for computing
 *           @f[
 *              \frac{
 *                  t_{ij} (t_{dk} \xi_k) - t_{dj} (t_{ik} \xi_k)
 *              }{
 *                  (t_{dk} \xi_k)^2
 *              }
 *           @f]
 *           where @f$ i, j \in \{0 \ldots d-1\} @f$ and @f$ k \in \{0 \ldots d\} @f$.
 */
template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransformDerivative : public ExpressionTraits<TValue> {
private:
    using TransformationMatrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1, Dimension+1>;

public:
    CIE_DEFINE_CLASS_POINTERS(ProjectiveTransformDerivative)

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

public:
    /// @brief Identity by default.
    ProjectiveTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(ConstSpan in, Span out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    static constexpr unsigned size() noexcept;

    /// @brief Compute the determinant of the projective transform's jacobian.
    TValue evaluateDeterminant(ConstSpan in) const;

private:
    friend class ProjectiveTransform<TValue,Dimension>;

    /// @brief Construct from a @ref ProjectiveTransform.
    ProjectiveTransformDerivative(Ref<const ProjectiveTransform<TValue,Dimension>> rProjection);

    /// @brief Construct directly from a transformation matrix representing a @ref ProjectiveTransform.
    ProjectiveTransformDerivative(Ref<const TransformationMatrix> rTransformationMatrix) noexcept;

private:
    TransformationMatrix _projectionMatrix;
}; // class ProjectiveTransformDerivative



/** @brief Class representing a projective transform.
 *  @details Uniquely defines a mapping between any pair of generalized quads.
 *  @note Unfortunately, I haven't found the proper term for quadrilaterals
 *        generalized to arbitrary dimensions, so @a generalized @a quad
 *        will have to do for now.
 */
template <concepts::Numeric TValue, unsigned Dimension>
class ProjectiveTransform : private ExpressionTraits<TValue> {
public:
    CIE_DEFINE_CLASS_POINTERS(ProjectiveTransform)

    using TransformationMatrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1, Dimension+1>;

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using Derivative = ProjectiveTransformDerivative<TValue,Dimension>;

    using Inverse = ProjectiveTransform;

    using Point = typename Kernel<Dimension,TValue>::Point;

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
     *  @param transformed Array of the transformed cube's vertices.
     */
    ProjectiveTransform(std::span<const Point> transformed);

    /// @brief Apply the transformation on a vector defined by the provided components.
    void evaluate(ConstSpan in, Span out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    static constexpr unsigned size() noexcept;

    /// @brief Construct the inverse transform.
    Inverse makeInverse() const;

    /// @brief Construct the derivative of the projective transform.
    Derivative makeDerivative() const;

    /// @brief Get the matrix representation of the transformation.
    Ref<const TransformationMatrix> getTransformationMatrix() const noexcept;

private:
    /// @brief Construct from a precomputed transformation matrix.
    ProjectiveTransform(RightRef<TransformationMatrix> rMatrix) noexcept;

    /// @brief Get the matrix representation of the transformation.
    Ref<TransformationMatrix> getTransformationMatrix() noexcept;

    /// @brief Compute the transformation matrix from the homogeneous representation of transformed points.
    /// @param[in] pTransformedBegin Ptr to the first component of the homogenized transformed points.
    /// @param[out] rMatrix Transformation matrix to write to.
    /// @warning @a pTransformedBegin is mutated during construction.
    static void computeTransformationMatrix(Ptr<TValue> pTransformedBegin,
                                            Ref<TransformationMatrix> rMatrix);

private:
    TransformationMatrix _transformationMatrix;
}; // class ProjectiveTransform


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/ProjectiveTransform_impl.hpp"
