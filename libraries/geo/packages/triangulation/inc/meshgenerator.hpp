#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <array>
#include <utility>
#include <numbers>
#include <span>


namespace cie::geo {

using Vector2D = std::array<double,2>;
using Vertex2D = Vector2D;
using Vertex2DVector = std::vector<Vertex2D>;

using IndexVector = std::vector<Size>;

using TriangleConnectivity = std::array<Size,3>;
using Triangulation = std::pair<Vertex2DVector, std::vector<TriangleConnectivity>>;

// Aggregates triangulation quality parameters
struct TriangulationParameters {
    double edgeLength    = 0.00;
    double goodChopRatio = 0.68;
    double divisionAngle = std::numbers::pi / 6.00;
}; // struct TriangulationParameters

Triangulation triangulate(
    std::span<const Vertex2D> polygon,
    TriangulationParameters parameters = {});

Triangulation triangulate(
    std::span<const Vertex2D> vertices,
    const std::vector<IndexVector>& polygonRegions,
    TriangulationParameters parameters = {});

} // namespace cie::geo
