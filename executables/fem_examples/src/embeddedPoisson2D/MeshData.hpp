#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/definitions.hpp"
#include "embeddedPoisson2D/CellData.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/MeshBase.hpp"
#include "packages/numeric/inc/KDTreeQuadraturePointFactory.hpp"
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"


namespace cie::fem {


/// @brief Data structure common to the entire @ref Graph "mesh".
class MeshData : public MeshBase<Ansatz> {
public:
    using DomainData = Scalar;

    MeshData();

    MeshData(
        RightRef<Ansatz> rAnsatzSpace,
        RightRef<std::vector<std::pair<DomainData,std::vector<Scalar>>>> domainTriangles,
        std::span<const std::pair<DomainData,Scalar>> domainMap,
        Ref<const cie::io::JSONObject> rConfiguration);

    KDTreeQuadraturePointFactory<
        Dimension,
        Scalar,
        Cell,
        Scalar
    > makeQuadratureRule(Ref<const Cell> rCell) const;

    void subdomain(
        std::span<const Scalar> points,
        std::span<DomainData> subdomains) const;

    std::span<const std::pair<
        DomainData,
        Scalar>
    > domainMap() const noexcept;

    void setCells(RightRef<std::vector<Cell>> rCells) noexcept;

    [[nodiscard]] constexpr std::span<const Cell> cells() const noexcept {
        return _cells;
    }

    [[nodiscard]] constexpr std::span<Cell> cells() noexcept {
        return _cells;
    }

private:
    friend struct io::GraphML::Serializer<MeshData>;

    friend struct io::GraphML::Deserializer<MeshData>;

    /// @brief Set of quadrature points for a default local hypercube.
    /// @details These quadrature points are used while constructing
    ///          cell-specific quadrature rules.
    std::vector<QuadraturePoint<Dimension,Scalar,Scalar>> _quadraturePointSet;

    std::vector<std::pair<
        DomainData,
        std::vector<Scalar>>
    > _domainTriangles;

    std::vector<std::pair<
        DomainData,
        Scalar>
    > _domainMap;

    std::vector<Cell> _cells;

    std::array<unsigned,2> _treeDepthRange;
}; // class MeshData


} // namespace cie::fem
