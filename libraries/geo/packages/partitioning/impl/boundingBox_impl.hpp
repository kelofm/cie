#pragma once

// help the language server
#include "packages/partitioning/inc/boundingBox.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::geo {


template <concepts::HasBoundingBox TObject>
inline const typename TObject::BoundingBox&
boundingBox(TObject& rObject) noexcept {
    return rObject.boundingBox();
}


template <unsigned Dimension, class TCoordinate>
stack::Box<Dimension,TCoordinate>
boundingBox(Ref<const stack::Box<Dimension,TCoordinate>> rObject) noexcept {
    return rObject;
}



template <concepts::StaticContainer TPoint>
requires std::is_same_v<TPoint,typename GetTraits<TPoint>::Type::Point>
inline const AABBox<GetTraits<TPoint>::Type::Dimension,typename GetTraits<TPoint>::Type::Coordinate>
boundingBox(const TPoint& rPoint) noexcept {
    constexpr unsigned Dimension = GetTraitsT<TPoint>::Dimension;
    using Coordinate = typename GetTraitsT<TPoint>::Coordinate;
    return AABBox<Dimension,Coordinate>(
        rPoint,
        detail::makeOrigin<Dimension,Coordinate>()
    );
}


} // namespace cie::geo
