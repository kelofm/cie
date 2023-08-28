#ifndef CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_HPP
#define CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_HPP

// --- Utilty Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class AffineTransform;


///@addtogroup fem
///@{


/// @brief Expression representing the derivative of @ref AffineTransform.
template <concepts::Numeric TValue, unsigned Dimension>
class AffineTransformDerivative : public ExpressionTraits<TValue>
{
private:
    using TransformationMatrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension,Dimension>;

public:
    CIE_DEFINE_CLASS_POINTERS(AffineTransformDerivative)

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity by default.
    AffineTransformDerivative() noexcept;

    /// @brief This function is constant so the arguments are not necessary.
    void evaluate(ConstIterator it_begin,
                  ConstIterator it_end,
                  Iterator it_out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Compute the determinant of the affine transform's jacobian.
    TValue evaluateDeterminant(ConstIterator it_begin, ConstIterator it_end) const;

private:
    friend class AffineTransform<TValue,Dimension>;

    AffineTransformDerivative(Ref<const AffineTransform<TValue,Dimension>> r_transform);

private:
    TransformationMatrix _matrix;
}; // class AffineTransformDerivative



/** @brief Class representing an affine transformation.
 *  @details Uniquely defines a mapping between any pair of triangles in the specified dimension.
 */
template <concepts::Numeric TValue, unsigned Dimension>
class AffineTransform final : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(AffineTransform)

    using TransformationMatrix = typename Kernel<Dimension,TValue>::dense::template static_matrix<Dimension+1, Dimension+1>;

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity transform by default
    AffineTransform() noexcept;

    /** @brief Affine transformation from (D+1) transformed points.
     *
     *  @details The transformation is characterized by how it deforms a cube
     *           defined on @f$ [-1,1]^D @f$. Not all corners of the deformed
     *           cube are necessary, but only the base @f$ [-1]^D @f$ and its
     *           adjecent ones (@f$ D+1 @f$ points in total). For example, the
     *           transformed location of the following vertices are required in 2D:
     *           @f[ \begin{bmatrix}
     *           -1 &  1 & -1 \\
     *           -1 & -1 &  1
     *           \end{bmatrix} @f]
     *
     *  @details To avoid unintended overlapping transformations, the order of vertices
     *           should match the following pattern (example in 3D):
     *           @f[ \begin{bmatrix}
     *           -1 &  1 & -1 & -1 \\
     *           -1 & -1 &  1 & -1 \\
     *           -1 & -1 & -1 &  1
     *           \end{bmatrix} @f]
     *
     *  @param it_transformedBegin iterator pointing to the transformed cube's base @f$ [-1]^D @f$.
     *  @param it_transformedEnd iterator past the last transformed point (should be identical to it_transformedBegin + D + 1 + 1).
     */
    template <concepts::Iterator TPointIt>
    AffineTransform(TPointIt it_transformedBegin,
                    TPointIt it_transformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components.
    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the affine transform.
    AffineTransformDerivative<TValue,Dimension> makeDerivative() const noexcept;

    /// @brief Construct the inverse transform.
    AffineTransform makeInverse() const;

    /// @brief Get the matrix representation of the transformation.
    Ref<const TransformationMatrix> getTransformationMatrix() const noexcept;

private:
    /// @brief Construct from a precomputed transformation matrix.
    AffineTransform(RightRef<TransformationMatrix> r_matrix) noexcept;

    /// @brief Get the matrix representation of the transformation.
    Ref<TransformationMatrix> getTransformationMatrix() noexcept;

    /// @brief Compute the transformation matrix from the homogeneous representation of transformed points.
    /// @param[in] p_transformedBegin Ptr to the first component of the homogenized transformed points.
    /// @param[out] r_matrix Transformation matrix to write to.
    static void computeTransformationMatrix(Ptr<const TValue> p_transformedBegin,
                                            Ref<TransformationMatrix> r_matrix);

private:
    TransformationMatrix _transformationMatrix;
};

///@}

} // namespace cie::fem::maths

#include "packages/maths/impl/AffineTransform_impl.hpp"

#endif
