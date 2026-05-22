#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"


namespace cie::fem {


/// @brief Data structure unique to @ref Graph::Edge "boundaries" between cells.
struct BoundaryData {
    static constexpr inline unsigned Dimension = 2;

    constexpr BoundaryData() noexcept = default;

    constexpr BoundaryData(BoundaryID boundary) noexcept
        : _boundary(boundary)
    {}

    /// @brief Boundary identifier of the shared boundary between the adjacent cells.
    constexpr BoundaryID boundary() const noexcept {
        return _boundary;
    }

private:
    BoundaryID _boundary;
}; // struct BoundaryData


} // namespace cie::fem
