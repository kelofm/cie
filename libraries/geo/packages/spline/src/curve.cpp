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

template <concepts::Numeric TValue, unsigned Dimension>
std::array<TValue,Dimension> deBoorOptimized(
    TValue parametricCoordinate,
    std::size_t iKnotSpan,
    std::size_t polynomialOrder,
    Ref<const std::span<const TValue>> knotVector,
    Ref<const std::array<std::span<const TValue>,Dimension>> coordinates,
    Ref<const std::array<std::span<TValue>,Dimension>> buffers) {
        std::array<TValue,Dimension> output;
        #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) CIE_CHECK(
            buffers[iDimension].size() == polynomialOrder + 1,
            std::format(
                "expecting buffer of size {} for coordinate component {}, but got {}",
                polynomialOrder + 1,
                iDimension,
                buffers[iDimension].size()))
        #endif

        for( std::size_t j = 0; j < polynomialOrder + 1; ++j ) {
            const std::size_t iPoint = j + iKnotSpan - polynomialOrder;
            for (unsigned iDim=0u; iDim<Dimension; ++iDim) {
                buffers[iDim][j] = coordinates[iDim][iPoint];
            } // for iDim in range(Dimension)
        } // for j in range(polynomialOrder + 1)

        for(std::size_t r = 1; r < polynomialOrder + 1; ++r) {
            for(std::size_t j = polynomialOrder; j > r - 1; --j) {
                const TValue tj = knotVector[j + iKnotSpan - polynomialOrder];
                const TValue alpha = (tj - parametricCoordinate) / (tj - knotVector[j + iKnotSpan + 1 - r]);
                const TValue beta = static_cast<TValue>(1) - alpha;
                for (const auto& rBuffer : buffers)
                    rBuffer[j] = alpha * rBuffer[j] + beta * rBuffer[j - 1];
            } // for j in range(polynomialOrder, r - 1, -1)
        } // for r in range(1, polynomialOrder + 1)

        std::transform(
            buffers.begin(),
            buffers.end(),
            output.begin(),
            [] (const auto& rSpan) {return rSpan.back();});
        return output;
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


template <concepts::Numeric TValue, unsigned Dimension>
void evaluateBSplineCurve(
    std::span<const TValue> tCoordinates,
    Ref<const std::array<std::span<const TValue>,Dimension>> controlPoints,
    std::span<const TValue> knotVector,
    Ref<const std::array<std::span<TValue>,Dimension>> output) {
        const std::size_t numberOfSamples = tCoordinates.size( );
        const std::size_t numberOfPoints = controlPoints.front().size( );
        const std::size_t m = knotVector.size( );
        const std::size_t polynomialOrder = m - numberOfPoints - 1;

        #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        for (unsigned iDim=0u; iDim<Dimension; ++iDim) {
            CIE_CHECK(
                controlPoints[iDim].size() == numberOfPoints,
                std::format(
                    "expecting coordinate input array for component {} to be of size {}, but got {}",
                    iDim,
                    controlPoints[iDim].size(),
                    numberOfPoints))
            CIE_CHECK(
                output[iDim].size() == numberOfSamples,
                std::format(
                    "expecting coordinate output array for component {} to be of size {}, but got {}",
                    iDim,
                    output[iDim].size(),
                    numberOfSamples))
        }
        #endif

        std::vector<TValue> deBoorBuffer(Dimension * (polynomialOrder + 1));
        std::array<std::span<TValue>,Dimension> bufferViews;
        for (unsigned iDim=0u; iDim<Dimension; ++iDim)
            bufferViews[iDim] = std::span<TValue>(
                deBoorBuffer.data() + iDim * (polynomialOrder + 1),
                polynomialOrder + 1);

        for(std::size_t i = 0; i < numberOfSamples; ++i) {
            const std::size_t iKnotSpan = findKnotSpan(
                tCoordinates[i],
                numberOfPoints,
                knotVector);
            const auto sample = deBoorOptimized<TValue,Dimension>(
                tCoordinates[i],
                iKnotSpan,
                polynomialOrder,
                knotVector,
                controlPoints,
                bufferViews);
            for (unsigned iDim=0u; iDim<Dimension; ++iDim)
                output[iDim][i] = sample[iDim];
        }
}


template void evaluateBSplineCurve<double,2>(
    std::span<const double>,
    Ref<const std::array<std::span<const double>,2>>,
    std::span<const double>,
    Ref<const std::array<std::span<double>,2>>);

template void evaluateBSplineCurve<float,2>(
    std::span<const float>,
    Ref<const std::array<std::span<const float>,2>>,
    std::span<const float>,
    Ref<const std::array<std::span<float>,2>>);

template void evaluateBSplineCurve<double,3>(
    std::span<const double>,
    Ref<const std::array<std::span<const double>,3>>,
    std::span<const double>,
    Ref<const std::array<std::span<double>,3>>);

template void evaluateBSplineCurve<float,3>(
    std::span<const float>,
    Ref<const std::array<std::span<const float>,3>>,
    std::span<const float>,
    Ref<const std::array<std::span<float>,3>>);


StaticArray<std::vector<double>, 2> evaluate2DCurveDeBoor(const std::vector<double>& tCoordinates,
                                                          const std::vector<double>& xCoordinates,
                                                          const std::vector<double>& yCoordinates,
                                                          const std::vector<double>& knotVector) {
    StaticArray<std::vector<double>,2> out {
        std::vector<double>(tCoordinates.size()),
        std::vector<double>(tCoordinates.size())};
    evaluateBSplineCurve<double,2>(
        tCoordinates,
        {xCoordinates, yCoordinates},
        knotVector,
        {out.front(), out.back()});
    return out;
}


StaticArray<std::vector<float>, 2> evaluate2DCurveDeBoor(const std::vector<float>& tCoordinates,
                                                         const std::vector<float>& xCoordinates,
                                                         const std::vector<float>& yCoordinates,
                                                         const std::vector<float>& knotVector) {
    StaticArray<std::vector<float>,2> out {
        std::vector<float>(tCoordinates.size()),
        std::vector<float>(tCoordinates.size())};
    evaluateBSplineCurve<float,2>(
        tCoordinates,
        {xCoordinates, yCoordinates},
        knotVector,
        {out.front(), out.back()});
    return out;
}

} // cie::geo