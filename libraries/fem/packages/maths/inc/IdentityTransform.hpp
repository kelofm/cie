#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class IdentityTransform : public ExpressionTraits<TValue> {
public:
    static constexpr inline unsigned ParametricDimension = Dimension;

    static constexpr inline unsigned PhysicalDimension = Dimension;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::BufferSpan;

    using Derivative = IdentityTransform;

    using Inverse = IdentityTransform;

public:
    IdentityTransform() noexcept = default;

    void evaluate(ConstSpan in, Span out, BufferSpan) const noexcept
    {std::copy(in.begin(), in.end(), out.begin());}

    static constexpr unsigned size() noexcept
    {return Dimension;}

    static constexpr unsigned bufferSize() noexcept
    {return 0u;}

    constexpr Inverse makeInverse() const noexcept
    {return *this;}

    constexpr Derivative makeDerivative() const noexcept
    {return *this;}

    constexpr TValue evaluateDeterminant([[maybe_unused]] ConstSpan in, BufferSpan) const noexcept
    {return static_cast<TValue>(1);}
}; // class IdentityTransform


} // namespace cie::fem::maths
