#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/maths/inc/power.hpp"

// --- GEO Includes ---
#include "packages/primitives/inc/Primitive.hpp"
#include "packages/primitives/inc/concepts.hpp"

// --- STL Includes ---
#include <tuple>
#include <span>


namespace cie::geo {


namespace stack {


template <unsigned Dimension, concepts::Numeric TCoordinate>
class Box : public Traits<Dimension,TCoordinate> {
public:
    constexpr Box() noexcept = default;

    constexpr Box(Ref<const std::span<const TCoordinate,Dimension>> rBase,
                  Ref<const std::span<const TCoordinate,Dimension>> rLengths) noexcept;

    constexpr Box(Ref<const StaticArray<TCoordinate,Dimension>> rBase,
                  Ref<const StaticArray<TCoordinate,Dimension>> rLengths) noexcept;

    constexpr Ref<const StaticArray<TCoordinate,Dimension>> base() const noexcept;

    constexpr std::span<TCoordinate,Dimension> base() noexcept;

    constexpr Ref<const StaticArray<TCoordinate,Dimension>> lengths() const noexcept;

    constexpr std::span<TCoordinate,Dimension> lengths() noexcept;

    constexpr bool at(Ref<const std::span<const TCoordinate,Dimension>> rPoint) const noexcept;

    constexpr bool at(Ref<const typename Box::Point> rPoint) const noexcept;

    template <unsigned iDim = 0u>
    static constexpr bool at(Ptr<const TCoordinate> pPointBegin,
                             Ptr<const TCoordinate> pBaseBegin,
                             Ptr<const TCoordinate> pLengthBegin) noexcept;

    static void makeCorners(Ref<const std::span<const TCoordinate,Dimension>> rBase,
                            Ref<const std::span<const TCoordinate,Dimension>> rLengths,
                            Ref<const std::span<TCoordinate,Dimension*intPow(2,Dimension)>> rCorners) noexcept;

private:
    StaticArray<TCoordinate,Dimension> _base, _lengths;
}; // class Box


} // namespace stack


/// @brief Box template.
template <Size Dimension, concepts::Numeric TCoordinate>
class Box : public AbsPrimitive<Dimension,TCoordinate>,
            public stack::Box<Dimension,TCoordinate>
{
private:
    using Base = AbsPrimitive<Dimension,TCoordinate>;

public:
    using typename Base::Point;

    using primitive_constructor_arguments = std::tuple<Point,Point>;

public:
    Box() noexcept;

    Box(const Point& rBase, const Point& rLengths) noexcept;

    template <class ContainerType1, class ContainerType2>
    requires concepts::Container<ContainerType1,TCoordinate>
             && concepts::Container<ContainerType2,TCoordinate>
    Box(const ContainerType1& rBase,
        const ContainerType2& rLengths) noexcept;

    virtual Bool isDegenerate() const override;
};


namespace boolean {


/// @brief Box with point membership test
template <Size Dimension, concepts::Numeric TCoordinate = Double>
class Box
    : public ::cie::geo::Box<Dimension,TCoordinate>,
      public Object<Dimension,Bool,TCoordinate>
{
public:
    Box() noexcept = default;

    template <class ContainerType1, class ContainerType2>
    requires concepts::Container<ContainerType1,TCoordinate>
             && concepts::Container<ContainerType2,TCoordinate>
    Box(const ContainerType1& rBase,
        const ContainerType2& rLengths) noexcept;

    virtual Bool at(const typename Box<Dimension,TCoordinate>::Point& rPoint) const override;
};


} // namespace boolean


} // namespace cie::geo


namespace cie::concepts {

template <class T>
concept Box
= Primitive<T>
  && requires ( T instance, const T constInstance )
{
    {instance.base()};
    {constInstance.base()};
    {instance.lengths()};
    {constInstance.lengths()};
};

} // namespace cie::concepts


#include "packages/primitives/impl/Box_impl.hpp"
