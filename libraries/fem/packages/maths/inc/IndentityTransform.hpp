#pragma once

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
    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

    using Derivative = IdentityTransform;

    using Inverse = IdentityTransform;

public:
    IdentityTransform() noexcept = default;

    void evaluate(ConstSpan in, Span out) const noexcept
    {std::copy(in.begin(), in.end(), out.begin());}

    constexpr unsigned size() const noexcept
    {return Dimension;}

    constexpr Inverse makeInverse() const noexcept
    {return *this;}

    constexpr Derivative makeDerivative() const noexcept
    {return *this;}

    constexpr TValue evaluateDeterminant([[maybe_unused]] ConstSpan in) const noexcept
    {return static_cast<TValue>(1);}
}; // class IdentityTransform


} // namespace cie::fem::maths
