#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <vector>


namespace cie::fem::maths {


template <class TValue, unsigned TParametricDimension, unsigned TPhysicalDimension>
class BSpline : public ExpressionTraits<TValue> {
public:
    constexpr inline static unsigned ParametricDimension = TParametricDimension;

    constexpr inline static unsigned PhysicalDimension = TPhysicalDimension;

    using typename ExpressionTraits<TValue>::Value;

    using typename ExpressionTraits<TValue>::Span;

    using typename ExpressionTraits<TValue>::ConstSpan;

    using typename ExpressionTraits<TValue>::BufferSpan;

    BSpline() noexcept = default;

    BSpline(
        std::span<const std::span<const TValue>,PhysicalDimension> controlPointCoordinates,
        std::span<const TValue> knots)
    requires (ParametricDimension == 1u);

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    static constexpr unsigned size() noexcept {
        return PhysicalDimension;
    }

    static constexpr unsigned bufferSize() noexcept {
        return 0u;
    }

private:
    constexpr std::size_t controlPointCount(unsigned iDimension) const noexcept {
        return _controlPointCount[iDimension];
    }

    std::span<const TValue> controlPointCoordinates(unsigned iDimension) const noexcept
    requires (ParametricDimension == 1u) {
        return {
            _data.data() + iDimension * this->controlPointCount(0),
            this->controlPointCount(0)};
    }

    std::span<const TValue> knots() const noexcept
    requires (ParametricDimension == 1u) {
        auto pBegin = _data.data() + PhysicalDimension * this->controlPointCount(0);
        return {
            pBegin,
            static_cast<std::size_t>(std::distance(pBegin, _data.data() + _data.size()))};
    }

    std::span<TValue> controlPointCoordinates(unsigned iDimension) noexcept
    requires (ParametricDimension == 1u)  {
        return {
            _data.data() + iDimension * this->controlPointCount(0),
            this->controlPointCount(0)};
    }

    std::span<TValue> knots() noexcept
    requires (ParametricDimension == 1u)  {
        auto pBegin = _data.data() + PhysicalDimension * this->controlPointCount(0);
        return {
            pBegin,
            static_cast<std::size_t>(std::distance(pBegin, _data.data() + _data.size()))};
    }

    /// @brief Vector containing all geometric properties.
    /// @details The stored data consists of PhysicalDimension + 1 arrays:
    ///          - control point x coordinates
    ///          - control point y coordinates
    ///          - ...
    ///          - knot vector
    std::vector<TValue> _data;

    std::array<std::size_t,ParametricDimension> _controlPointCount;
}; // class BSpline


} // namespace cie::fem::maths
