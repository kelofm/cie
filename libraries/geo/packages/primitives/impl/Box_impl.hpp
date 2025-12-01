#pragma once

// help the language server
#include "packages/primitives/inc/Box.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/maths/inc/OuterProduct.hpp" // OuterProduct

// --- STL Includes ---
#include <algorithm>
#include <limits>


namespace cie::geo {


namespace stack {


template <unsigned Dimension, concepts::Numeric TCoordinate>
constexpr Box<Dimension,TCoordinate>::Box(Ref<const std::span<const TCoordinate,Dimension>> rBase,
                                          Ref<const std::span<const TCoordinate,Dimension>> rLengths) noexcept
    : _base(),
      _lengths()
{
    std::copy_n(rBase.data(), Dimension, _base.data());
    std::copy_n(rLengths.data(), Dimension, _lengths.data());
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
constexpr Box<Dimension,TCoordinate>::Box(Ref<const StaticArray<TCoordinate,Dimension>> rBase,
                                          Ref<const StaticArray<TCoordinate,Dimension>> rLengths) noexcept
    : _base(rBase),
      _lengths(rLengths)
{}


template <unsigned Dimension, concepts::Numeric TCoordinate>
StaticArray<TCoordinate,Dimension>
constexpr Box<Dimension,TCoordinate>::base() const noexcept {
    return _base;
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
std::span<TCoordinate,Dimension>
constexpr Box<Dimension,TCoordinate>::base() noexcept {
    return std::span<TCoordinate,Dimension>(_base.data(), Dimension);
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
StaticArray<TCoordinate,Dimension>
constexpr Box<Dimension,TCoordinate>::lengths() const noexcept {
    return _lengths;
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
std::span<TCoordinate,Dimension>
constexpr Box<Dimension,TCoordinate>::lengths() noexcept {
    return std::span<TCoordinate,Dimension>(_lengths.data(), Dimension);
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
StaticArray<TCoordinate,Dimension>
constexpr Box<Dimension,TCoordinate>::opposite() const noexcept {
    StaticArray<TCoordinate,Dimension> out;
    std::transform(
        _base.begin(),
        _base.end(),
        _lengths.begin(),
        out.begin(),
        std::plus<TCoordinate>());
    return out;
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
bool
constexpr Box<Dimension,TCoordinate>::at(Ref<const std::span<const TCoordinate,Dimension>> rPoint) const noexcept {
    return this->at(rPoint.data(), _base.data(), _lengths.data());
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
bool
constexpr Box<Dimension,TCoordinate>::at(Ref<const typename Box::Point> rPoint) const noexcept {
    return this->at(rPoint.data(), _base.data(), _lengths.data());
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
template <unsigned iDim>
constexpr bool
Box<Dimension,TCoordinate>::at(Ptr<const TCoordinate> pPointBegin,
                               Ptr<const TCoordinate> pBaseBegin,
                               Ptr<const TCoordinate> pLengthBegin) noexcept {
    if constexpr (iDim != Dimension) {
        const bool lessThanLowerBound = *pPointBegin < *pBaseBegin;
        const bool lessThanUpperBound = *pPointBegin < (*pBaseBegin) + (*pLengthBegin);
        if (lessThanLowerBound == lessThanUpperBound) return false;

        if constexpr (iDim + 1 < Dimension) {
            if (!Box::template at<iDim+1>(++pPointBegin, ++pBaseBegin, ++pLengthBegin)) return false;
        }
    }
    return true;
}


template <unsigned Dimension, concepts::Numeric TCoordinate>
void
Box<Dimension,TCoordinate>::makeCorners(Ref<const std::span<const TCoordinate,Dimension>> rBase,
                                        Ref<const std::span<const TCoordinate,Dimension>> rLengths,
                                        Ref<const std::span<TCoordinate,Dimension*intPow(2,Dimension)>> rCorners) noexcept {
    StaticArray<std::uint8_t,Dimension> state;
    std::fill_n(state.data(), Dimension, static_cast<std::uint8_t>(0));
    std::uint8_t iCorner = 0;

    do {
        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
            Ref<TCoordinate> rCornerComponent = rCorners[iCorner * Dimension + iDimension];
            rCornerComponent = rBase[iDimension];
            if (state[iDimension]) rCornerComponent += rLengths[iDimension];
        } // for iDimension in range(Dimension)

        ++iCorner;
    } while (::cie::maths::OuterProduct<Dimension>::next(2u, state.data()));
}


} // namespace stack


/* --- Box --- */

template < Size Dimension,
           concepts::Numeric TCoordinate >
Box<Dimension,TCoordinate>::Box(const Point& rBase, const Point& rLengths) noexcept
{
    #ifdef CIE_ENABLE_RUNTIME_GEOMETRY_CHECKS
    for (const auto& length : rLengths)
        CIE_DEBUG_CHECK(0 <= length, "Edge lengths of a box must be non-negative (" << length << ")")
    #endif
    std::copy_n(rBase.data(), Dimension, this->base().data());
    std::copy_n(rLengths.data(), Dimension, this->lengths().data());
}


template < Size Dimension,
           concepts::Numeric TCoordinate >
template <class ContainerType1, class ContainerType2>
requires concepts::Container<ContainerType1,TCoordinate>
         && concepts::Container<ContainerType2,TCoordinate>
Box<Dimension,TCoordinate>::Box(const ContainerType1& rBase,
                                const ContainerType2& rLengths) noexcept
{
    CIE_OUT_OF_RANGE_CHECK( rBase.size() == Dimension )
    CIE_OUT_OF_RANGE_CHECK( rLengths.size() == Dimension )

    #ifdef CIE_ENABLE_RUNTIME_GEOMETRY_CHECKS
    for (const auto& length : rLengths)
        CIE_DEBUG_CHECK(0 <= length, "Edge lengths of a box must be non-negative")
    #endif

    std::copy(rBase.begin(),
              rBase.end(),
              this->base().begin() );
    std::copy(rLengths.begin(),
              rLengths.end(),
              this->lengths().begin() );
}


template < Size Dimension,
           concepts::Numeric TCoordinate >
Box<Dimension,TCoordinate>::Box() noexcept
    : Box<Dimension,TCoordinate>(detail::makeOrigin<Dimension,TCoordinate>(),
                                 detail::makeOrigin<Dimension,TCoordinate>())
{
}


template < Size Dimension,
           concepts::Numeric TCoordinate >
Bool
Box<Dimension,TCoordinate>::isDegenerate() const
{
    Bool degenerate = false;
    for (const auto component : this->lengths())
        if (component < std::numeric_limits<TCoordinate>::epsilon()) {
            degenerate = true;
            break;
        }

    return degenerate;
}


/* --- boolean::Box --- */

namespace boolean {


template <  Size Dimension,
            concepts::Numeric TCoordinate   >
template <class ContainerType1, class ContainerType2>
    requires concepts::Container<ContainerType1,TCoordinate>
                && concepts::Container<ContainerType2,TCoordinate>
Box<Dimension,TCoordinate>::Box(const ContainerType1& rBase,
                                const ContainerType2& rLengths) noexcept
    : cie::geo::Box<Dimension,TCoordinate>(rBase, rLengths)
{
}


template <Size Dimension, concepts::Numeric TCoordinate>
Bool Box<Dimension,TCoordinate>::at(const typename Box<Dimension,TCoordinate>::Point& rPoint) const
{
    const auto& rBase = this->base();
    const auto& rLengths = this->lengths();
    return stack::Box<Dimension,TCoordinate>::at(
        rPoint.data(),
        rBase.data(),
        rLengths.data());
}


} // namespace boolean


} // namespace cie::geo
