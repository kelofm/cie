// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- Internal Includes ---
#include "packages/spline/inc/basisfunctions.hpp"
#include "packages/spline/inc/curve.hpp"

// --- STL Includes ---
#include <algorithm>
#include <cmath>
#include <string>
#include <format>


namespace cie::geo {


template <concepts::Numeric TValue>
std::size_t findKnotSpan(TValue t,
                         std::size_t numberOfControlPoints,
                         Ref<const std::span<const TValue>> knotVector) {
    TValue tolerance = 1e-10;

    // Check if t resides within the allowed bounds
    if(t < knotVector.front() || knotVector.back() < t) {
        CIE_THROW(
            OutOfRangeException,
            "t out range: t = " + std::to_string( t ) +
            " but can only be within " + std::to_string( knotVector.front( ) ) +
            " and " + std::to_string( knotVector.back( ) ) + "\n"
        )
    }

    if( std::abs( t - knotVector[numberOfControlPoints + 1] ) < tolerance ) {
        return numberOfControlPoints - 1;
    }

    auto result = std::upper_bound( knotVector.begin( ), knotVector.end( ), t );
    return std::distance(knotVector.begin( ), result - 1);
}


StaticArray<std::vector<double>, 2> evaluate2DCurve(const std::vector<double>& tCoordinates,
                                                    const std::vector<double>& xCoordinates,
                                                    const std::vector<double>& yCoordinates,
                                                    const std::vector<double>& knotVector) {
    std::size_t numberOfSamples = tCoordinates.size( );
    std::size_t numberOfPoints = xCoordinates.size( );
    std::size_t m = knotVector.size( );
    std::size_t p = m - numberOfPoints - 1;

    if( yCoordinates.size( ) != numberOfPoints )
        CIE_THROW( OutOfRangeException, "Inconsistent size in evaluate2DCurve" )

    std::vector<double> curveX( numberOfSamples, 0.0 );
    std::vector<double> curveY( numberOfSamples, 0.0 );

    for( std::size_t i = 0; i < numberOfSamples; ++i ) {
        double t = tCoordinates[i];

        for( std::size_t j = 0; j < numberOfPoints; ++j ) {
            double N = evaluateBSplineBasis( t, j, p, knotVector );
            curveX[i] += N * xCoordinates[j];
            curveY[i] += N * yCoordinates[j];
        }
    }

    return { curveX, curveY };
}

template <concepts::Numeric TValue>
StaticArray<TValue, 2> deBoorOptimized(TValue t,
                                       std::size_t knotSpanIndex,
                                       std::size_t polynomialDegree,
                                       Ref<const std::span<const TValue>> knotVector,
                                       Ref<const std::span<const TValue>> xCoordinates,
                                       Ref<const std::span<const TValue>> yCoordinates,
                                       std::span<TValue> dx,
                                       std::span<TValue> dy) {
    CIE_CHECK(
        dx.size() == polynomialDegree + 1,
        std::format(
            "expecting buffer of size {}, but got {}",
            polynomialDegree + 1,
            dx.size()))
    CIE_CHECK(
        dy.size() == polynomialDegree + 1,
        std::format(
            "expecting buffer of size {}, but got {}",
            polynomialDegree + 1,
            dy.size()))

    for( std::size_t j = 0; j < polynomialDegree + 1; ++j ) {
        dx[j] = xCoordinates[j + knotSpanIndex - polynomialDegree];
        dy[j] = yCoordinates[j + knotSpanIndex - polynomialDegree];
    }

    for(std::size_t r = 1; r < polynomialDegree + 1; ++r) {
        for(std::size_t j = polynomialDegree; j > r - 1; --j) {
            TValue tj = knotVector[j + knotSpanIndex - polynomialDegree];
            TValue alpha = ( tj - t ) / ( tj - knotVector[j + knotSpanIndex + 1 - r] );
            dx[j] = ( 1.0 - alpha ) * dx[j - 1] + alpha * dx[j];
            dy[j] = ( 1.0 - alpha ) * dy[j - 1] + alpha * dy[j];
        }
    }

    return { dx[polynomialDegree], dy[polynomialDegree] };
}

StaticArray<double, 2> deBoor( double t,
                              std::size_t knotSpanIndex,
                              std::size_t polynomialDegree,
                              const std::vector<double>& knotVector,
                              const std::vector<double>& xCoordinates,
                              const std::vector<double>& yCoordinates,
                              std::size_t refinementLevel )
{
    if( refinementLevel == polynomialDegree + 1 )
    {
        return { xCoordinates[knotSpanIndex], yCoordinates[knotSpanIndex] };
    }

    double a = ( t - knotVector[knotSpanIndex] ) / ( knotVector[knotSpanIndex + refinementLevel] - knotVector[knotSpanIndex] );

    StaticArray<double, 2> P1 = deBoor( t, knotSpanIndex - 1, polynomialDegree, knotVector, xCoordinates, yCoordinates, refinementLevel + 1 );
    StaticArray<double, 2> P2 = deBoor( t, knotSpanIndex, polynomialDegree, knotVector, xCoordinates, yCoordinates, refinementLevel + 1 );

    double Px = ( 1.0 - a ) * P1[0] + a * P2[0];
    double Py = ( 1.0 - a ) * P1[1] + a * P2[1];

    return { Px, Py};
}


template <concepts::Numeric TValue>
void evaluate2DCurveDeBoor(std::span<const TValue> tCoordinates,
                           std::span<const TValue> xCoordinates,
                           std::span<const TValue> yCoordinates,
                           std::span<const TValue> knotVector,
                           std::span<TValue> xOut,
                           std::span<TValue> yOut) {
    const std::size_t numberOfSamples = tCoordinates.size( );
    const std::size_t numberOfPoints = xCoordinates.size( );
    const std::size_t m = knotVector.size( );
    const std::size_t polynomialOrder = m - numberOfPoints - 1;

    CIE_CHECK(
        xOut.size() == numberOfSamples,
        "x coordinate output size " << xOut.size() << " does not match curve parameter array size " << numberOfSamples)
    CIE_CHECK(
        yOut.size() == numberOfSamples,
        "y coordinate output size " << yOut.size() << " does not match curve parameter array size " << numberOfSamples)

    if( yCoordinates.size( ) != numberOfPoints )
        CIE_THROW( OutOfRangeException, "Inconsistent size in evaluate2DCurveDeBoor" );

    std::vector<TValue> deBoorBuffer(2 * (polynomialOrder + 1));
    const std::span<TValue>
        xBuffer(
            deBoorBuffer.data(),
            polynomialOrder + 1),
        yBuffer(
            deBoorBuffer.data() + polynomialOrder + 1,
            polynomialOrder + 1);

    for(std::size_t i = 0; i < numberOfSamples; ++i) {
        const std::size_t iKnotSpan = findKnotSpan(
            tCoordinates[i],
            numberOfPoints,
            knotVector);
        const StaticArray<TValue, 2> Point = deBoorOptimized(
            tCoordinates[i],
            iKnotSpan,
            polynomialOrder,
            knotVector,
            xCoordinates,
            yCoordinates,
            xBuffer,
            yBuffer);
        xOut[i] = Point[0];
        yOut[i] = Point[1];
    }
}


StaticArray<std::vector<double>, 2> evaluate2DCurveDeBoor(const std::vector<double>& tCoordinates,
                                                          const std::vector<double>& xCoordinates,
                                                          const std::vector<double>& yCoordinates,
                                                          const std::vector<double>& knotVector) {
    StaticArray<std::vector<double>,2> out {
        std::vector<double>(tCoordinates.size()),
        std::vector<double>(tCoordinates.size())};
    evaluate2DCurveDeBoor<double>(
        tCoordinates,
        xCoordinates,
        yCoordinates,
        knotVector,
        out.front(),
        out.back());
    return out;
}


StaticArray<std::vector<float>, 2> evaluate2DCurveDeBoor(const std::vector<float>& tCoordinates,
                                                         const std::vector<float>& xCoordinates,
                                                         const std::vector<float>& yCoordinates,
                                                         const std::vector<float>& knotVector) {
    StaticArray<std::vector<float>,2> out {
        std::vector<float>(tCoordinates.size()),
        std::vector<float>(tCoordinates.size())};
    evaluate2DCurveDeBoor<float>(
        tCoordinates,
        xCoordinates,
        yCoordinates,
        knotVector,
        out.front(),
        out.back());
    return out;
}

} // cie::geo