#ifndef CIE_FEM_SCALE_TRANSLATE_TRANSFORM_HPP
#define CIE_FEM_SCALE_TRANSLATE_TRANSFORM_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransform;


/// @addtogroup fem
/// @{


template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransformDerivative : public ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(ScaleTranslateTransformDerivative)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity by default.
    ScaleTranslateTransformDerivative() noexcept;

    /// @brief Evaluate the derivative at the provided point.
    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const noexcept;

    /// @brief Evaluate the determinant of the original transform's jacobian at the provided location.
    TValue evaluateDeterminant(ConstIterator it_argumentBegin,
                               ConstIterator it_argumentEnd) const noexcept;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

private:
    friend class ScaleTranslateTransform<TValue,Dimension>;

    ScaleTranslateTransformDerivative(Ref<const ScaleTranslateTransform<TValue,Dimension>> r_transform) noexcept;

private:
    StaticArray<TValue,Dimension> _scales;
}; // class ScaleTranslateTransformDerivative



/// @brief Class representing independent scaling along orthogonal coordinate axes, followed by a translation.
/// @details Uniquely defines a mapping between axis-aligned hyperrectangles in D-dimensional space.
/// @note In order to keep an efficient representation and evalutation, this transformation does not provide
///       a means to invert it, since that would involve flipping the order of scaling and translating.
template <concepts::Numeric TValue, unsigned Dimension>
class ScaleTranslateTransform : private ExpressionTraits<TValue>
{
public:
    CIE_DEFINE_CLASS_POINTERS(ScaleTranslateTransform)

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

public:
    /// @brief Identity transform by default.
    ScaleTranslateTransform() noexcept;

    /** @brief Construct from the transformed base and the vertex opposite @f$ [1]^D @f$.
     *  @details The coordinates of the input transformed vertex are identical to the
     *           scaling coefficients of the transform, after undoing the translation.
     *  @note The number of input components must match the dimension.
     */
    template <concepts::Iterator TPointIt>
    ScaleTranslateTransform(TPointIt it_transformedBegin,
                            TPointIt it_transformedEnd);

    /// @brief Apply the transformation on a vector defined by the provided components
    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const;

    /// @brief Get the number of components written by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the transform.
    ScaleTranslateTransformDerivative<TValue,Dimension> makeDerivative() const noexcept;

private:
    friend class ScaleTranslateTransformDerivative<TValue,Dimension>;

    StaticArray<TValue,Dimension> _scales;

    StaticArray<TValue,Dimension> _offset;
}; // class ScaleTranslateTransform


/// @}


} // namespace cie::fem::maths

#include "packages/maths/impl/ScaleTranslateTransform_impl.hpp"

#endif
