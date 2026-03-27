#pragma once

// --- External Includes ---
#include "tsl/robin_set.h"
#include "tsl/robin_map.h"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>
#include <iterator>


namespace std {


inline bool operator==(pair<cie::fem::BoundaryID,cie::Size> left,
                       pair<cie::fem::BoundaryID,cie::Size> right) noexcept {
    return (left.first == right.first) && (left.second == right.second);
}


inline bool operator!=(pair<cie::fem::BoundaryID,cie::Size> left,
                       pair<cie::fem::BoundaryID,cie::Size> right) noexcept {
    return left.first != right.first && left.second != right.second;
}


template <>
struct hash<pair<cie::fem::BoundaryID,cie::Size>> {
    auto operator()(pair<cie::fem::BoundaryID,cie::Size> item) const
    {
        const auto tmp = hash<cie::fem::BoundaryID>()(item.first);
        return tmp ^ (hash<cie::Size>()(item.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); ///< from boost::hash_combine
    }
};


} // namespace std


namespace cie::fem {


template <maths::Expression TAnsatzSpace, cie::concepts::CallableWith<BoundaryID,Size> TFunctor>
void scanConnectivities(Ref<const TAnsatzSpace> rAnsatzSpace,
                        TFunctor&& rFunctor,
                        Ptr<const typename TAnsatzSpace::Value> pSampleBegin,
                        Ptr<const typename TAnsatzSpace::Value> pSampleEnd,
                        typename TAnsatzSpace::Value tolerance) {
    CIE_BEGIN_EXCEPTION_TRACING
    using Value = typename TAnsatzSpace::Value;
    constexpr unsigned Dimension = TAnsatzSpace::Dimension;
    static_assert(0 < Dimension);
    std::vector<typename TAnsatzSpace::Value> ansatzBuffer(rAnsatzSpace.bufferSize());

    const Size ansatzSize = rAnsatzSpace.size();
    const Size numberOfSamples = std::distance(pSampleBegin, pSampleEnd);

    StaticArray<Value,Dimension> argument;
    StaticArray<unsigned,Dimension-1> argumentState;
    DynamicArray<Value> valueBuffer(ansatzSize);

    std::fill_n(
        argumentState.data(),
        Dimension - 1,
        0u);

    constexpr unsigned maxBoundaries = 2 * Dimension;
    BoundaryID boundaryID;
    tsl::robin_set<std::pair<BoundaryID,Size>> activeBoundaries;

    for (unsigned iBoundary=0; iBoundary<maxBoundaries; ++iBoundary, ++boundaryID) {
        activeBoundaries.clear();
        const unsigned iDim = boundaryID.getDimension();

        do {
            // Compute the current sample point
            Size iState = 0;
            for (unsigned i=0; i<Dimension; ++i) {
                if (i == iDim) {
                    argument[i] = static_cast<Value>(boundaryID.getDirection() ? 1 : -1);
                } else {
                    argument[i] = *(pSampleBegin + argumentState[iState++]);
                }
            } // for i in range(Dimension)

            // Evaluate the ansatz space
            rAnsatzSpace.evaluate(
                argument,
                valueBuffer,
                ansatzBuffer);

            for (Size iAnsatz=0; iAnsatz<ansatzSize; ++iAnsatz) {
                if (tolerance < valueBuffer[iAnsatz]) {
                    activeBoundaries.emplace(boundaryID, iAnsatz);
                }
            }
        } while (cie::maths::OuterProduct<Dimension-1>::next(numberOfSamples, argumentState.data()));

        for (auto pair : activeBoundaries) {
            rFunctor(pair.first, pair.second);
        }
    } // for iBoundary in range(maxBoundaries)
    CIE_END_EXCEPTION_TRACING
}


template <unsigned Dimension>
template <maths::Expression TAnsatzSpace>
AnsatzMap<Dimension>::AnsatzMap(
    Ref<const TAnsatzSpace> rAnsatzSpace,
    std::size_t integrationOrder,
    utils::Comparison<typename TAnsatzSpace::Value> comparison)
    : _ansatzCount(rAnsatzSpace.size()),
      _topology() {
        CIE_BEGIN_EXCEPTION_TRACING

        using Value = typename TAnsatzSpace::Value;
        GaussLegendreQuadrature<Value> quadrature(integrationOrder, comparison);
        std::vector<typename TAnsatzSpace::Value> ansatzBuffer(rAnsatzSpace.bufferSize());

        if (quadrature.nodes().size() && _ansatzCount) {
            // Loop through viable boundary pairs.
            // These consist of oriented boundaries whose axes coincide with
            // the original local axes but may have flipped directions. Boundary
            // pairs' axes must point in the same direction in all directions but
            // their normals, where they must point in opposite directions.
            StaticArray<StaticArray<BoundaryID,2>,Dimension> axisOptions;
            StaticArray<Ptr<const BoundaryID>,Dimension> axisOptionBegins, axisOptionEnds;
            for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                axisOptions[iDimension][0] = BoundaryID(iDimension, false);
                axisOptions[iDimension][1] = BoundaryID(iDimension, true);
                axisOptionBegins[iDimension] = axisOptions[iDimension].data();
                axisOptionEnds[iDimension] = axisOptions[iDimension].data() + axisOptions[iDimension].size();
            } // for iDimension in range(Dimension)

            StaticArray<BoundaryID,Dimension> axisState;
            std::transform(
                axisOptions.begin(),
                axisOptions.end(),
                axisState.begin(),
                [](const auto& rArray){return rArray.front();});

            do {
                // Container for sample points on the boundary.
                // Since all but one axes are aligned, the locations of
                // coincident sample points are identical in both systems,
                // so we don't have to generate separate points and evaluate
                // them for the two boundaries. This also means that the only
                // coincident pairs of ansatz functions on each boundary are
                // the ones that don't vanish there.
                DynamicArray<StaticArray<Value,Dimension>> samplePoints;

                for (unsigned iBoundaryAxis=0u; iBoundaryAxis<Dimension; ++iBoundaryAxis) {
                    samplePoints.clear();

                    // Construct entries for the two boundaries along the current axis
                    // in the internal connectivity map.
                    typename ConnectivityMap::iterator
                        itNegativeBoundaryConnectivities,
                        itPositiveBoundaryConnectivities;
                    {
                        typename ConnectivityMap::mapped_type emptyConnectivities;

                        OrientedBoundary<Dimension> negativeBoundary, positiveBoundary;
                        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                            negativeBoundary[iDimension] = axisState[iDimension];
                            positiveBoundary[iDimension] = axisState[iDimension];
                        }
                        negativeBoundary.id() = BoundaryID(iBoundaryAxis, false);
                        positiveBoundary.id() = BoundaryID(iBoundaryAxis, true);

                        OrientedBoundary<Dimension> negativeLeft = negativeBoundary,
                                                    negativeRight = negativeBoundary;
                        negativeRight[iBoundaryAxis] = -BoundaryID(negativeRight[iBoundaryAxis]);
                        if (!_topology.emplace(
                                std::make_pair(negativeLeft, negativeRight),
                                emptyConnectivities
                            ).second) {
                            continue;
                        }

                        OrientedBoundary<Dimension> positiveLeft = positiveBoundary,
                                                    positiveRight = positiveBoundary;
                        positiveRight[iBoundaryAxis] = -BoundaryID(positiveRight[iBoundaryAxis]);
                        itPositiveBoundaryConnectivities = _topology.emplace(
                            std::make_pair(positiveLeft, positiveRight),
                            emptyConnectivities
                        ).first;

                        itNegativeBoundaryConnectivities = _topology.find(std::make_pair(negativeLeft, negativeRight));
                        CIE_OUT_OF_RANGE_CHECK(itNegativeBoundaryConnectivities != _topology.end())
                    }

                    // N-d index of the sample point.
                    StaticArray<Size,Dimension-1> argumentState;
                    std::fill_n(argumentState.data(), Dimension - 1, 0ul);

                    // Loop through the sample points as one-hot outer product of the
                    // provided sample coordinates. Coordinates of the boundary axis are
                    // omitted, since they may only take a value of -1 (negative boundary)
                    // or 1 (positive boundary).
                    do {
                        samplePoints.emplace_back();
                        unsigned iState = 0u;
                        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                            if (iDimension != iBoundaryAxis) {
                                samplePoints.back()[iDimension] = quadrature.nodes()[argumentState[iState++]];
                            }
                        }
                    } while (cie::maths::OuterProduct<Dimension-1>::next(quadrature.numberOfNodes(), argumentState.data()));

                    // Indices of ansatz functions that don't vanish on the current boundary.
                    // Begin by assuming that all functions vanish.
                    DynamicArray<bool> positiveSideVanish(_ansatzCount, true), negativeSideVanish(_ansatzCount, true);
                    tsl::robin_map<Size,tsl::robin_set<Size>> positiveSideCoincidents, negativeSideCoincidents;
                    positiveSideCoincidents.reserve(_ansatzCount);
                    negativeSideCoincidents.reserve(_ansatzCount);
                    {
                        tsl::robin_set<Size> all;
                        all.reserve(_ansatzCount);
                        for (Size iAnsatz=0ul; iAnsatz<_ansatzCount; ++iAnsatz) {
                            all.insert(iAnsatz);
                        }
                        for (Size iAnsatz=0ul; iAnsatz<_ansatzCount; ++iAnsatz) {
                            positiveSideCoincidents.emplace(iAnsatz, all);
                            negativeSideCoincidents.emplace(iAnsatz, all);
                        }
                    }

                    // Buffer to store the basis functions evaluated at a sample point.
                    DynamicArray<Value> valueBuffer(_ansatzCount);

                    // Buffer to store the integral of the squared difference of the ansatz functions.
                    DynamicArray<Value> squaredDifferenceIntegral(_ansatzCount);
                    std::fill_n(
                        squaredDifferenceIntegral.begin(),
                        _ansatzCount,
                        static_cast<Value>(0));

                    // Loop through all sample points and evaluate all basis functions at
                    // them. Then check which ones coincide and erase the pairs from the
                    // map that don't.
                    for (auto& rSamplePoint : samplePoints) {
                        // For practical reasons, we're not looping over boundaries, but
                        // axes, which means we have to check two boundaries:
                        // - one on the positive side of the axis
                        // - and another one on the opposite side.
                        // So the component of the current sample point related to the current
                        // axis must be set to 1 or -1, depending on which boundary we're checking.

                        // Negative side
                        rSamplePoint[iBoundaryAxis] = -1;
                        rAnsatzSpace.evaluate(rSamplePoint, valueBuffer, ansatzBuffer);
                        for (Size iLeftAnsatz=0u; iLeftAnsatz<_ansatzCount; ++iLeftAnsatz) {
                            negativeSideVanish[iLeftAnsatz] = negativeSideVanish[iLeftAnsatz] && comparison.equal(valueBuffer[iLeftAnsatz], 0);

                            for (Size iRightAnsatz=iLeftAnsatz; iRightAnsatz<_ansatzCount; ++iRightAnsatz) {
                                if (!comparison.equal(valueBuffer[iLeftAnsatz], valueBuffer[iRightAnsatz])) {
                                    negativeSideCoincidents[iLeftAnsatz].erase(iRightAnsatz);
                                    negativeSideCoincidents[iRightAnsatz].erase(iLeftAnsatz);
                                } // if left != right
                            } // for iRightAnsatz in range(iLeftAnsatz, _ansatzCount)
                        } // for iLeftAnsatz in range(_ansatzCount)

                        // Positive side
                        rSamplePoint[iBoundaryAxis] = 1;
                        rAnsatzSpace.evaluate(rSamplePoint, valueBuffer, ansatzBuffer);
                        for (Size iLeftAnsatz=0u; iLeftAnsatz<_ansatzCount; ++iLeftAnsatz) {
                            positiveSideVanish[iLeftAnsatz] = positiveSideVanish[iLeftAnsatz] && comparison.equal(valueBuffer[iLeftAnsatz], 0);

                            for (Size iRightAnsatz=iLeftAnsatz; iRightAnsatz<_ansatzCount; ++iRightAnsatz) {
                                if (!comparison.equal(valueBuffer[iLeftAnsatz], valueBuffer[iRightAnsatz])) {
                                    positiveSideCoincidents[iLeftAnsatz].erase(iRightAnsatz);
                                    positiveSideCoincidents[iRightAnsatz].erase(iLeftAnsatz);
                                } // if left != right
                            } // for iRightAnsatz in range(iLeftAnsatz, _ansatzCount)
                        } // for iLeftAnsatz in range(_ansatzCount)
                    } // for rSamplePoint in samplePoints

                    // Check which functions haven't vanished and are coincident,
                    // and add them to the internal map.
                    for (unsigned iAnsatz=0u; iAnsatz<_ansatzCount; ++iAnsatz) {
                        if (!negativeSideVanish[iAnsatz]) {
                            const auto it = negativeSideCoincidents.find(iAnsatz);
                            CIE_OUT_OF_RANGE_CHECK(it != negativeSideCoincidents.end())

                            DynamicArray<std::pair<Size,Size>> value;
                            value.reserve(it->second.size());
                            std::transform(
                                it->second.begin(),
                                it->second.end(),
                                std::back_inserter(value),
                                [iAnsatz](const Size iOtherAnsatz) -> std::pair<Size,Size> {
                                    return std::make_pair(iAnsatz, iOtherAnsatz);
                                });

                            std::copy(
                                value.begin(),
                                value.end(),
                                std::back_inserter(itNegativeBoundaryConnectivities.value()));
                        } // if !negativeSideVanish[iAnsatz]

                        if (!positiveSideVanish[iAnsatz]) {
                            const auto it = positiveSideCoincidents.find(iAnsatz);
                            CIE_OUT_OF_RANGE_CHECK(it != positiveSideCoincidents.end())

                            DynamicArray<std::pair<Size,Size>> value;
                            value.reserve(it->second.size());
                            std::transform(it->second.begin(),
                                        it->second.end(),
                                        std::back_inserter(value),
                                        [iAnsatz](const Size iOtherAnsatz) -> std::pair<Size,Size> {
                                                return std::make_pair(iAnsatz, iOtherAnsatz);
                                        });

                            std::copy(value.begin(),
                                    value.end(),
                                    std::back_inserter(itPositiveBoundaryConnectivities.value()));
                        } // if !positiveSideVanish[iAnsatz]
                    } // for iAnsatz in range(_ansatzCount)
                } // for iBoundaryAxis in range(dimension)
            } while (cie::maths::OuterProduct<Dimension>::next(
                axisOptionBegins.data(),
                axisOptionEnds.data(),
                axisState.data()));
        } // if sampleSize && _ansatzCount
        CIE_END_EXCEPTION_TRACING
}


template <maths::Expression TAnsatzSpace, class TValue>
requires std::is_same_v<TValue,typename TAnsatzSpace::Value>
AnsatzMap<TAnsatzSpace::Dimension>
makeAnsatzMap(Ref<const TAnsatzSpace> rAnsatzSpace,
              std::size_t integrationOrder,
              utils::Comparison<TValue> comparison) {
    return AnsatzMap<TAnsatzSpace::Dimension>(rAnsatzSpace, integrationOrder, comparison);
}


} // namespace cie::fem
