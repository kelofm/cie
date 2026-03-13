#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <vector>
#include <span>
#include <array>

namespace cie::geo {

/** @brief Evaluate B-Spline curve by summing up basis functions times control points.
 *  @param tCoordinates Parametric coordinates at which the curve is evaluated.
 *  @param xCoordinates The @p x coordinates of the control points.
 *  @param yCoordinates The @p y coordinates of the control points.
 *  @param knotVector Parametric coordinates of the control points.
 *  @return A vector of @p x and a vector of @p y coordinates with one value for each parametric
 *          coordinate tCoordinates.
 */
StaticArray<std::vector<double>, 2> evaluate2DCurve(
    const std::vector<double>& tCoordinates,
    const std::vector<double>& xCoordinates,
    const std::vector<double>& yCoordinates,
    const std::vector<double>& knotVector);

template <concepts::Numeric TValue, unsigned Dimension>
void evaluateBSplineCurve(
    std::span<const TValue> tCoordinates,
    Ref<const std::array<std::span<const TValue>,Dimension>> controlPoints,
    std::span<const TValue> knotVector,
    Ref<const std::array<std::span<TValue>,Dimension>> output);

//! Identical to evaluate2DCurve, but using De Boor's algorithm.
StaticArray<std::vector<double>, 2> evaluate2DCurveDeBoor(
    const std::vector<double>& tCoordinates,
    const std::vector<double>& xCoordinates,
    const std::vector<double>& yCoordinates,
    const std::vector<double>& knotVector );

//! Identical to evaluate2DCurve, but using De Boor's algorithm.
StaticArray<std::vector<float>, 2> evaluate2DCurveDeBoor(
    const std::vector<float>& tCoordinates,
    const std::vector<float>& xCoordinates,
    const std::vector<float>& yCoordinates,
    const std::vector<float>& knotVector );

/*! De Boor's algorithm for evaluating (x, y) at one parametric coordinate t. The parameter
 *  recursionLevel has a default value of 1, which will be used if no argument is passed. */
StaticArray<double, 2> deBoor(
    double t,
    std::size_t knotSpanIndex,
    std::size_t polynomialDegree,
    const std::vector<double>& knotVector,
    const std::vector<double>& xCoordinates,
    const std::vector<double>& yCoordinates,
    std::size_t recursionLevel = 1 );

//! Same as above but without recursion.
template <concepts::Numeric TValue, unsigned Dimension>
std::array<TValue,Dimension> deBoorOptimized(
    TValue parametricCoordinate,
    std::size_t iKnotSpan,
    std::size_t polynomialOrder,
    Ref<const std::span<const TValue>> knotVector,
    Ref<const std::array<std::span<const TValue>,Dimension>> coordinates,
    Ref<const std::array<std::span<TValue>,Dimension>> buffers);

//! Determines the knot span of the parametric coordinate t.
template <concepts::Numeric TValue>
std::size_t findKnotSpan(
    TValue t,
    std::size_t numberOfControlPoints,
    Ref<const std::span<const TValue>> knotVector);

} // cie::geo

