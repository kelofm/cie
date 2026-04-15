#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/numeric/inc/CellBase.hpp"
#include "packages/numeric/inc/MeshBase.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/BoxBoundable.hpp"
#include "packages/primitives/inc/Object.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"


namespace cie::fem {


/// @brief Data structure unique to each @ref Graph::Vertex "cell".
struct CellData
    : public geo::BoxBoundable<Dimension,Scalar>,
      public CellBase<Dimension,Scalar,SpatialTransform,Scalar> {
    using BoxBase = geo::BoxBoundable<cie::fem::Dimension,Scalar>;

    using CellBase = CellBase<cie::fem::Dimension,Scalar,SpatialTransform,Scalar>;

    CellData() noexcept = default;

    CellData(
        VertexID id,
        AnsatzID ansatzID,
        Scalar diffusivity,
        OrientedAxes<Dimension> axes,
        RightRef<SpatialTransform> rSpatialTransform) noexcept;

    Scalar diffusivity() const noexcept {
        return this->data();
    }

    bool at(geo::BoxBoundable<Dimension,Scalar>::Point point) const;

protected:
    void computeBoundingBoxImpl(BoundingBox& rBox) noexcept override;
}; // struct CellData


static_assert(::cie::concepts::SamplableGeometry<CellData>);


template <>
struct io::GraphML::Serializer<CellData>
    : io::GraphML::Serializer<CellData::CellBase>
{};


template <>
struct io::GraphML::Deserializer<CellData>
    : io::GraphML::Deserializer<CellData::CellBase>
{};


} // namespace cie::fem
