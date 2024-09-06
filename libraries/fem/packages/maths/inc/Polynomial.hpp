#ifndef CIE_FEM_MATHS_POLYNOMIAL_HPP
#define CIE_FEM_MATHS_POLYNOMIAL_HPP

// --- FEM Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::fem::maths {


/// @brief @ref Expression representing a scalar polynomial.
template <class TValue>
class Polynomial : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::Iterator;

    using typename ExpressionTraits<TValue>::ConstIterator;

    using Derivative = Polynomial;

    using Coefficients = DynamicArray<TValue>;

public:
    /// @brief Uninitialized by default.
    Polynomial() noexcept = default;

    /// @brief Construct from a container of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    Polynomial(RightRef<Coefficients> rCoefficients) noexcept;

    /// @brief Construct from a range of coefficients.
    /// @details The input coefficients are expected to be sorted
    ///          in the order of their corresponding monomials.
    template <concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
    Polynomial(TItBegin itBegin, TItEnd itEnd);

    Polynomial(Polynomial&& rRhs) noexcept = default;

    Polynomial(Ref<const Polynomial> rRhs) = default;

    Ref<Polynomial> operator=(RightRef<Polynomial> rRhs) noexcept = default;

    Ref<Polynomial> operator=(Ref<const Polynomial> rRhs) = default;

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itResultBegin) const;

    /// @brief Get the number of scalar components returned by @ref evaluate.
    unsigned size() const noexcept;

    /// @brief Construct the derivative of the @ref Polynomial.
    /// @details The returned object is also a @ref Polynomial,
    ///          irrespective of the current polynomial order.
    Polynomial makeDerivative() const;

protected:
    Coefficients _coefficients;
}; // class Polynomial


} // namespace cie::fem::maths

#include "packages/maths/impl/Polynomial_impl.hpp"

#endif
