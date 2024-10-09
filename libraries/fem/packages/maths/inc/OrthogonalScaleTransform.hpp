#ifndef CIE_FEM_ORTHOGONAL_SCALE_TRANSFORM_HPP
#define CIE_FEM_ORTHOGONAL_SCALE_TRANSFORM_HPP

// --- FEM Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransform;


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransformDerivative : public ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(OrthogonalScaleTransformDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity by default.
    OrthogonalScaleTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const noexcept;

    /// @brief Evaluate the determinant of the original transform's jacobian at the provided location.
    TValue evaluateDeterminant(ConstIterator itArgumentBegin,
                               ConstIterator itArgumentEnd) const noexcept;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

private:
    friend class OrthogonalScaleTransform<TValue,Dimension>;

    OrthogonalScaleTransformDerivative(Ref<const OrthogonalScaleTransform<TValue,Dimension>> rTransform) noexcept;

private:
    StaticArray<TValue,Dimension> _scales;
}; // class OrthogonalScaleTransformDerivative



/// @brief Class representing independent scaling along orthogonal coordinate axes.
/// @details Uniquely defines a mapping between axis-aligned hyperrectangles in
///          D-dimensional space, that have their base vertices at the origin.
template <concepts::Numeric TValue, unsigned Dimension>
class OrthogonalScaleTransform : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(OrthogonalScaleTransform)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

    using Derivative = OrthogonalScaleTransformDerivative<TValue,Dimension>;

    using Inverse = OrthogonalScaleTransform;

public:
    /// @brief Identity transform by default.
    OrthogonalScaleTransform() noexcept;

    /** @brief Construct from the transformed vertex opposite the base @f$ [1]^D @f$.
     *  @details The coordinates of the input transformed vertex are identical to the
     *           scaling coefficients of the transform.
     *  @note The number of input components must match the dimension.
     */
    template <concepts::Iterator TPointIt>
    OrthogonalScaleTransform(TPointIt itTransformedBegin,
                             TPointIt itTransformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components
    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the transform.
    Derivative makeDerivative() const noexcept;

    /// @brief Construct the inverse transform.
    Inverse makeInverse() const;

private:
    friend class OrthogonalScaleTransformDerivative<TValue,Dimension>;

    StaticArray<TValue,Dimension> _scales;
}; // class OrthogonalScaleTransform


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/OrthogonalScaleTransform_impl.hpp"

#endif
