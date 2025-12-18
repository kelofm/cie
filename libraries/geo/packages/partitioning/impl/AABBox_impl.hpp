#pragma once

// help the language server
#include "packages/partitioning/inc/AABBox.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/OuterProduct.hpp"


namespace cie::geo {


template <Size Dimension, concepts::Numeric TCoordinate>
Bool AABBox<Dimension,TCoordinate>::contains(const stack::Box<Dimension,TCoordinate>& rBox) const noexcept
{
    // All points inside
    for (Size dim=0; dim<Dimension; ++dim) {
        if (rBox.base()[dim] < this->base()[dim]
            ||
            this->base()[dim] + this->lengths()[dim] < rBox.base()[dim] + rBox.lengths()[dim])
            return false;
    }

    return true;
}


template <Size Dimension, concepts::Numeric TCoordinate>
Bool AABBox<Dimension,TCoordinate>::intersects(const stack::Box<Dimension,TCoordinate>& rBox) const noexcept
{
    //if (rBox.contains(*this)) return true;

    bool hasCornerOutside = false;
    bool hasCornerInside = false;

    // At least one point outside and one inside
    StaticArray<std::uint8_t,Dimension> samplePointState;
    std::fill_n(samplePointState.data(), Dimension, static_cast<std::uint8_t>(0));
    do {
        typename AABBox::Point corner;
        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
            // Construct corner.
            corner[iDimension] = rBox.base()[iDimension];
            if (samplePointState[iDimension]) {
                corner[iDimension] += rBox.lengths()[iDimension];
            }
        } // for iDimension in range(Dimension)

        // Check whether the corner is inside or outside the bbox.
        if (this->at(corner)) hasCornerInside = true;
        else hasCornerOutside = true;

        if (hasCornerInside && hasCornerOutside) return true;
    } while (cie::maths::OuterProduct<Dimension>::next(2u, samplePointState.data()));

    return false;
}


template <Size Dimension, concepts::Numeric TCoordinate>
void AABBox<Dimension,TCoordinate>::include(const stack::Box<Dimension,TCoordinate>& rBox) noexcept
{
    for (Size dim=0; dim<Dimension; ++dim) {
        auto dBase = rBox.base()[dim] - this->base()[dim];

        if (dBase < static_cast<TCoordinate>(0)) {
            this->base()[dim]    += dBase;
            this->lengths()[dim] -= dBase;
        }

        auto thisMax = this->base()[dim] + this->lengths()[dim];
        auto boxMax  = rBox.base()[dim] + rBox.lengths()[dim];

        if (thisMax < boxMax)
            this->lengths()[dim] = boxMax - this->base()[dim];
    } // for dim in range(Dimension)
}


} // namespace cie::geo
