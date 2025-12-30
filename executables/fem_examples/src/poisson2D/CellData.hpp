#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/numeric/inc/Cell.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/BoxBoundable.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"


namespace cie::fem {


/// @brief Data structure unique to each @ref Graph::Vertex "cell".
struct CellData
    : public geo::BoxBoundable<Dimension,Scalar>,
      public CellBase<Dimension,Scalar,SpatialTransform,Scalar> {
    using BoxBase = geo::BoxBoundable<cie::fem::Dimension,Scalar>;

    using CellBase = CellBase<cie::fem::Dimension,Scalar,SpatialTransform,Scalar>;

    using CellBase::Dimension;

    using LocalCoordinate = Kernel<Dimension,Scalar>::LocalCoordinate;

    using GlobalCoordinate = Kernel<Dimension,Scalar>::GlobalCoordinate;

    CellData() noexcept = default;

    CellData(VertexID id,
             unsigned short iAnsatz,
             Scalar diffusivity,
             OrientedAxes<Dimension> axes,
             RightRef<SpatialTransform> rSpatialTransform) noexcept
        : BoxBase(),
          CellBase(
            id,
            iAnsatz,
            axes,
            std::move(rSpatialTransform),
            std::move(diffusivity))
    {}

    Scalar diffusivity() const noexcept {
        return this->data();
    }

    bool at(geo::BoxBoundable<Dimension,Scalar>::Point point) const {
        const utils::Comparison<Scalar> comparison(1e-8, 1e-10);
        StaticArray<LocalCoordinate,Dimension> local;
        this->transform(
            Kernel<Dimension,Scalar>::cast<GlobalCoordinate>(std::span<const Scalar,Dimension>(point.data(), Dimension)),
            Kernel<Dimension,Scalar>::view(local));

        return std::all_of(
            local.begin(),
            local.end(),
            [&comparison](Scalar coordinate) {
                return comparison.less(std::abs(coordinate), static_cast<Scalar>(1))
                    || comparison.equal(std::abs(coordinate), static_cast<Scalar>(1));}
        );
    }

protected:
    void computeBoundingBoxImpl(BoundingBox& rBox) noexcept override {
        BoundingBox::Point opposite;
        StaticArray<LocalCoordinate,Dimension> localCorner;
        StaticArray<GlobalCoordinate,Dimension> globalCorner;

        std::fill_n(
            rBox.base().data(),
            Dimension,
            std::numeric_limits<Scalar>::max());
        std::fill_n(
            opposite.data(),
            Dimension,
            std::numeric_limits<Scalar>::lowest());

        StaticArray<std::uint8_t,Dimension> state {0u, 0u};
        StaticArray<LocalCoordinate,2> ordinates {-1.0, 1.0};

        do {
            // Compute the corner in local space.
            std::transform(state.begin(),
                           state.end(),
                           localCorner.begin(),
                           [&ordinates](std::uint8_t iOrdinate) {
                                return ordinates[iOrdinate];
                           });

            // Transform the corner to global space.
            this->transform(
                Kernel<Dimension,Scalar>::view(localCorner),
                Kernel<Dimension,Scalar>::view(globalCorner));

            // Extend box definition.
            for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                rBox.base()[iDimension] = std::min<Scalar>(rBox.base()[iDimension], globalCorner[iDimension]);
                opposite[iDimension] = std::max<Scalar>(opposite[iDimension], globalCorner[iDimension]);
            } // for iDimension in range(Dimension)
        } while (cie::maths::OuterProduct<Dimension>::next(2u, state.data()));

        std::transform(opposite.begin(),
                       opposite.end(),
                       rBox.base().begin(),
                       rBox.lengths().begin(),
                       std::minus<Scalar>());
    }
}; // struct CellData


template <>
struct io::GraphML::Serializer<CellData>
    : io::GraphML::Serializer<CellData::CellBase>
{};


template <>
struct io::GraphML::Deserializer<CellData>
    : io::GraphML::Deserializer<CellData::CellBase>
{};


} // namespace cie::fem
