#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/GraphTraits.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/numeric/inc/CellBase.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/BoxBoundable.hpp"
#include "packages/primitives/inc/Object.hpp"


namespace cie::fem {


/// @details Vertices in the mesh's graph do not store additional information.
///          Cells are referenced by their corresponding vertices' IDs in the
///          graph of the mesh.
using CellData = void;


/// @brief Data structure unique to each @ref Graph::Vertex "cell".
class Cell
    :   public geo::BoxBoundable<Dimension,Scalar>,
        public CellBase<Dimension,Scalar,SpatialTransform,Scalar> {
public:
    using BoxBase = geo::BoxBoundable<cie::fem::Dimension,Scalar>;

    using CellBase = CellBase<cie::fem::Dimension,Scalar,SpatialTransform,Scalar>;

    Cell() noexcept = default;

    Cell(
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
}; // struct Cell


static_assert(::cie::concepts::SamplableGeometry<Cell>);


} // namespace cie::fem
