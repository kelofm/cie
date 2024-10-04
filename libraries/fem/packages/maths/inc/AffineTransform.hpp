#ifndef CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_HPP
#define CIE_FEM_MATHS_ABS_AFFINE_TRANSFORM_HPP

// --- Utilty Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"

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
    void evaluate(ConstIterator itBegin,
                  ConstIterator itEnd,
                  Iterator itOut) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Compute the determinant of the affine transform's jacobian.
    TValue evaluateDeterminant(ConstIterator itBegin, ConstIterator itEnd) const;

private:
    friend class AffineTransform<TValue,Dimension>;

    AffineTransformDerivative(Ref<const AffineTransform<TValue,Dimension>> rTransform);

private:
    TransformationMatrix _matrix;
}; // class AffineTransformDerivative



/** @brief Class representing an affine transformation.
 *  @details Uniquely defines a mapping between any pair of triangles in the specified dimension.
 *           Implements the @ref SpatialTransform interface.
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

    using Derivative = AffineTransformDerivative<TValue,Dimension>;

    using Inverse = AffineTransform;

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
     *  @param itTransformedBegin iterator pointing to the transformed cube's base @f$ [-1]^D @f$.
     *  @param itTransformedEnd iterator past the last transformed point (should be identical to itTransformedBegin + D + 1 + 1).
     */
    template <concepts::Iterator TPointIt>
    AffineTransform(TPointIt itTransformedBegin,
                    TPointIt itTransformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components.
    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the affine transform.
    Derivative makeDerivative() const noexcept;

    /// @brief Construct the inverse transform.
    Inverse makeInverse() const;

    /// @brief Get the matrix representation of the transformation.
    Ref<const TransformationMatrix> getTransformationMatrix() const noexcept;

private:
    /// @brief Construct from a precomputed transformation matrix.
    AffineTransform(RightRef<TransformationMatrix> rMatrix) noexcept;

    /// @brief Get the matrix representation of the transformation.
    Ref<TransformationMatrix> getTransformationMatrix() noexcept;

    /// @brief Compute the transformation matrix from the homogeneous representation of transformed points.
    /// @param[in] pTransformedBegin Ptr to the first component of the homogenized transformed points.
    /// @param[out] rMatrix Transformation matrix to write to.
    static void computeTransformationMatrix(Ptr<const TValue> pTransformedBegin,
                                            Ref<TransformationMatrix> rMatrix);

private:
    TransformationMatrix _transformationMatrix;
};

///@}

} // namespace cie::fem::maths

#include "packages/maths/impl/AffineTransform_impl.hpp"

#endif
