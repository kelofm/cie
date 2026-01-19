#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"

// --- Utility Includes ---
//#include "packages/concurrency/inc/sycl.hpp"


namespace cie::fem {


/// @brief Radius of the circle Dirichlet conditions are imposed on.
constexpr double boundaryRadius             = 2.5e-1;

/// @brief Minimum depth of the spatial tree used to find intersections between domain cells and boundary cells.
constexpr unsigned minBoundaryTreeDepth     = 3;

/// @brief Maximum depth of the spatial tree used to find intersections between domain cells and boundary cells.
constexpr unsigned maxBoundaryTreeDepth     = 20;

/// @brief Minimum norm of a boundary segment to integrate over.
constexpr double minBoundarySegmentNorm     = 1e-12;

/// @brief Number of spatial dimensions the problem is defined on.
constexpr unsigned Dimension = 2u;

/// @brief Floating point scalar type to use.
using Scalar = double;

/// @brief Spatial transform type mapping cells' local space to global space.
/// @details The flexibility of this transform directly defines what kind
///          of geometries the cells can represent in a boundary conforming manner.
using SpatialTransform = maths::ScaleTranslateTransform<Scalar,Dimension>;

/// @brief Cells' basis functions' type.
using Basis = maths::PolynomialView<Scalar>;

/// @brief Cells' ansatz space type.
/// @details Spans the full outer product space of the basis functions.
using Ansatz = maths::AnsatzSpaceView<Basis,Dimension>;
using AnsatzDerivative = maths::AnsatzSpaceDerivativeView<Basis,Dimension>;


} // namespace cie::fem
