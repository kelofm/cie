#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/io/inc/Traits.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
class IdentityTransformDerivative
    :   public ExpressionTraits<TValue>,
        public cie::io::TriviallySerializableBase {
public:
    static constexpr inline unsigned ParametricDimension = Dimension;

    static constexpr inline unsigned PhysicalDimension = Dimension;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::BufferSpan;

    using Inverse = IdentityTransformDerivative;

    template <template <class ...> class, class>
    using Rebind = IdentityTransformDerivative;

public:
    constexpr IdentityTransformDerivative() noexcept = default;

    constexpr void evaluate(ConstSpan, Span out, BufferSpan) const noexcept {
        auto it = out.begin();
        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension)
            for (unsigned jDimension=0u; jDimension<Dimension; ++jDimension)
                *it++ = iDimension == jDimension;
    }

    static constexpr unsigned size() noexcept
    {return Dimension * Dimension;}

    static constexpr unsigned bufferSize() noexcept
    {return 0u;}

    constexpr Inverse makeInverse() const noexcept
    {return *this;}

    constexpr TValue evaluateDeterminant(ConstSpan, BufferSpan) const noexcept
    {return static_cast<TValue>(1);}

    void serialize(
        Ref<cie::io::Traits::SerializerStream>,
        tags::Binary = {}) const
    {}

    template <class TAllocator>
    void deserialize(
        Ref<cie::io::Traits::DeserializerStream>,
        Ref<IdentityTransformDerivative>,
        Ref<const TAllocator>,
        tags::Binary = {})
    {}
}; // class IdentityTransformDerivative


template <concepts::Numeric TValue, unsigned Dimension>
class IdentityTransform
    :   public ExpressionTraits<TValue>,
        public cie::io::TriviallySerializableBase {
public:
    static constexpr inline unsigned ParametricDimension = Dimension;

    static constexpr inline unsigned PhysicalDimension = Dimension;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::BufferSpan;

    using Derivative = IdentityTransformDerivative<TValue,Dimension>;

    using Inverse = IdentityTransform;

    template <template <class ...> class, class>
    using Rebind = IdentityTransform;

public:
    constexpr IdentityTransform() noexcept = default;

    constexpr void evaluate(ConstSpan in, Span out, BufferSpan) const noexcept
    {std::copy(in.begin(), in.end(), out.begin());}

    static constexpr unsigned size() noexcept
    {return Dimension;}

    static constexpr unsigned bufferSize() noexcept
    {return 0u;}

    constexpr Derivative makeDerivative() const noexcept
    {return Derivative();}

    constexpr Inverse makeInverse() const noexcept
    {return Inverse();}

    void serialize(
        Ref<cie::io::Traits::SerializerStream>,
        tags::Binary = {}) const
    {}

    template <class TAllocator>
    void deserialize(
        Ref<cie::io::Traits::DeserializerStream>,
        Ref<IdentityTransform>,
        Ref<const TAllocator>,
        tags::Binary = {})
    {}
}; // class IdentityTransform


} // namespace cie::fem::maths
