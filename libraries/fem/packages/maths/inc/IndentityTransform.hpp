#ifndef CIE_FEM_MATHS_IDENTITY_TRANSFORM_HPP
#define CIE_FEM_MATHS_IDENTITY_TRANSFORM_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class IdentityTransform : public ExpressionTraits<TValue>
{
public:
    using typename ExpressionTraits<TValue>::ConstIterator;

    using typename ExpressionTraits<TValue>::Iterator;

    using Derivative = IdentityTransform;

public:
    IdentityTransform() noexcept = default;

    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const noexcept
    {std::fill(it_out, it_out + Dimension, static_cast<TValue>(1));}

    constexpr unsigned size() const noexcept
    {return Dimension;}

    constexpr IdentityTransform makeInverse() const noexcept
    {return *this;}

    constexpr Derivative makeDerivative() const noexcept
    {return *this;}
}; // class IdentityTransform


} // namespace cie::fem::maths


#endif
