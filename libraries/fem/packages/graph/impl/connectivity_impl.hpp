#ifndef CIE_FEM_GRAPH_CONNECTIVITY_IMPL_HPP
#define CIE_FEM_GRAPH_CONNECTIVITY_IMPL_HPP

// --- External Includes ---
#include "tsl/robin_set.h"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/maths/inc/CartesianProduct.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <algorithm>
#include <iterator>


namespace std {


inline bool operator==(pair<cie::fem::BoundaryID,cie::Size> left,
                       pair<cie::fem::BoundaryID,cie::Size> right) noexcept
{
    return (left.first == right.first) && (left.second == right.second);
}


inline bool operator!=(pair<cie::fem::BoundaryID,cie::Size> left,
                       pair<cie::fem::BoundaryID,cie::Size> right) noexcept
{
    return left.first != right.first && left.second != right.second;
}


template <>
struct hash<pair<cie::fem::BoundaryID,cie::Size>>
{
    auto operator()(pair<cie::fem::BoundaryID,cie::Size> item) const
    {
        const auto tmp = hash<cie::fem::BoundaryID>()(item.first);
        return tmp ^ (hash<cie::Size>()(item.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
    }
};


} // namespace std


namespace cie::fem {


template <maths::Expression TAnsatzSpace, concepts::CallableWith<BoundaryID,Size> TFunctor>
void scanConnectivities(Ref<const TAnsatzSpace> rAnsatzSpace,
                        TFunctor&& rFunctor,
                        Ptr<const typename TAnsatzSpace::Value> pSampleBegin,
                        Ptr<const typename TAnsatzSpace::Value> pSampleEnd,
                        typename TAnsatzSpace::Value tolerance)
{
    CIE_BEGIN_EXCEPTION_TRACING
    using Value = typename TAnsatzSpace::Value;
    constexpr unsigned Dimension = TAnsatzSpace::Dimension;
    static_assert(0 < Dimension);

    const Size ansatzSize = rAnsatzSpace.size();
    const Size numberOfSamples = std::distance(pSampleBegin, pSampleEnd);

    StaticArray<Value,Dimension> argument;
    StaticArray<Size,Dimension-1> argumentState;
    DynamicArray<Value> valueBuffer(ansatzSize);

    std::fill(argumentState.begin(),
              argumentState.end(),
              static_cast<Size>(0));

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
            rAnsatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data());

            for (Size iAnsatz=0; iAnsatz<ansatzSize; ++iAnsatz) {
                if (tolerance < valueBuffer[iAnsatz]) {
                    activeBoundaries.emplace(boundaryID, iAnsatz);
                }
            }
        } while (maths::CartesianProduct<Dimension-1>::next(numberOfSamples, argumentState.data()));

        for (auto pair : activeBoundaries) {
            rFunctor(pair.first, pair.second);
        }
    } // for iBoundary in range(maxBoundaries)
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
template <maths::Expression TAnsatzSpace>
requires (std::is_same_v<typename TAnsatzSpace::Value,TValue>)
AnsatzMap<TValue>::AnsatzMap(Ref<const TAnsatzSpace> rAnsatzSpace,
                             std::span<const TValue> samples,
                             utils::Comparison<TValue> comparison)
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr unsigned Dimension = TAnsatzSpace::Dimension;
    static_assert(0 < Dimension);

    using Value = typename TAnsatzSpace::Value;
    const Size ansatzSize = rAnsatzSpace.size();   // <== total number of ansatz functions
    const Size numberOfSamples = samples.size();    // <== number of sample nodes per dimension

    StaticArray<Value,Dimension> argument; // <== sample point
    StaticArray<Size,Dimension-1> argumentState; // <== sample point index state for the outer product

    // Buffer to evaluate the ansatz functions to
    // - the first half [0;ansatzSize] stores values on the negative boundary
    // - the second half [ansatzSize;2*ansatzSize] stores values on the positive boundary
    DynamicArray<Value> valueBuffer(2 * ansatzSize);

    // Reset the sample point index state before iterating over the outer product
    std::fill(argumentState.begin(),
              argumentState.end(),
              static_cast<unsigned>(0));

    tsl::robin_map<
        Size,                   // <== ansatz function index on the negative boundary
        tsl::robin_set<Size>    // <== indices of ansatz functions with identical values on the positive boundary
    > ansatzPairs;

    // Indices of ansatz functions that don't vanish on the given boundary
    tsl::robin_set<std::pair<BoundaryID,Size>> nonzeros;

    // Loop over dimensions and find coincident ansatz functions on opposite boundaries
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        // Assume every ansatz function vanishes on both boundaries
        nonzeros.clear();

        // Assume every ansatz function is coincident with every other function
        // (consistent with the assumption of all functions vanishing on both boundaries)
        {
            tsl::robin_set<Size> all;
            all.reserve(ansatzSize);
            for (Size i=0; i<ansatzSize; ++i) all.insert(i);
            for (Size iAnsatz=0; iAnsatz<ansatzSize; ++iAnsatz) ansatzPairs.insert_or_assign(iAnsatz, all);
        }

        do {
            // Compute the current sample point
            unsigned iState = 0;
            for (unsigned i=0; i<Dimension; ++i) {
                if (i != iDim) {
                    argument[i] = samples[argumentState[iState++]];
                }
            } // for i in range(Dimension)

            // Evaluate the ansatz space on the negative boundary
            argument[iDim] = Value(-1);
            rAnsatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data());

            // Evaluate the ansatz space on the positive boundary
            argument[iDim] = Value(1);
            rAnsatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data() + ansatzSize);

            // Compare ansatz function values at opposite boundaries
            for (unsigned iAnsatz=0; iAnsatz<ansatzSize; ++iAnsatz) {
                // Check whether the negative boundary is nonzero
                if (!comparison.equal(valueBuffer[iAnsatz], Value(0))) {
                    nonzeros.emplace(BoundaryID(iDim, false), iAnsatz);
                }

                // Check whether the positive boundary is nonzero
                if (!comparison.equal(valueBuffer[iAnsatz + ansatzSize], Value(0))) {
                    nonzeros.emplace(BoundaryID(iDim, true), iAnsatz);
                }
            }

            // Update the coincidence map
            for (auto it=ansatzPairs.begin(); it!=ansatzPairs.end(); ++it) {
                const Size iNegative = it.key(); // <== index of the ansatz function on the negative boundary
                DynamicArray<Size> erase;

                for (Size iPositive : it.value()) {
                    if (!comparison.equal(valueBuffer[iNegative], valueBuffer[iPositive + ansatzSize])) {
                        erase.push_back(iPositive);
                    } // if not coincident
                } // for iPositive in pair.second

                for (auto iErase : erase) {
                    it.value().erase(iErase);
                }
            } // for pair in ansatzPairs
        } while (maths::CartesianProduct<Dimension-1>::next(numberOfSamples, argumentState.data()));

        // Finished looping over every sample point on both boundaries
        // in the current direction, and evaluating all ansatz functions
        // on them. What we ended up with is:
        // - ansatzPairs: a map associating each ansatz function index with
        //                a set of other ansatz function indices coincident
        //                on the opposite boundaries.
        // - nonzeros: an array storing which ansatz functions don't vanish
        //             on what boundary.

        // Erase entries from the coincidence map that vanish on the negative boundary
        for (Size iAnsatz=0; iAnsatz<ansatzSize; ++iAnsatz) {
            if (nonzeros.find({{iDim, false}, iAnsatz}) == nonzeros.end()) {
                ansatzPairs.erase(iAnsatz);
            }
        }

        // Check whether each ansatz function that doesn't vanish on the negative
        // boundary has exactly ONE coincident ansatz function on the positive boundary.
        for (auto it=ansatzPairs.begin(); it!=ansatzPairs.end(); ++it) {
            if (it.value().size() != 1) {
                std::stringstream message;
                message << "Ansatz function " << it.key() << " at the negative boundary of dimension " << iDim
                        << " should have exactly 1 coincident ansatz function on the positive boundary, but has "
                        << it.value().size();
                if (!it.value().empty()) {
                    auto p=it.value().begin();
                    message << "(" << *p++;
                    for (; p!=it.value().end(); ++p) {
                        message << ", " << *p;
                    }
                    message << ")";
                }
                message << ".\n";
                CIE_THROW(Exception, message.str())
            }
        }

        // Check whether each ansatz funtion that doesn't vanish on the positive boundary
        // has a coincident ansatz function on the negative boundary.
        for (auto [boundary, iAnsatz] : nonzeros) {
            if (boundary.getDirection()) {
                bool foundCoincidentPair = false;
                for (const auto& rPair : ansatzPairs) {
                    if (rPair.second.find(iAnsatz) != rPair.second.end()) {
                        foundCoincidentPair = true;
                        break;
                    }
                }
                if (!foundCoincidentPair) {
                    std::stringstream message;
                    message << "Ansatz function " << iAnsatz
                            << " does not vanish on the positive boundary of dimension " << iDim
                            << " but has no coincident ansatz function on the negative boundary.\n";
                    CIE_THROW(Exception, message.str())
                }
            }
        }

        // Update the connectivity map
        DynamicArray<std::pair<Size,Size>> connectivities;
        connectivities.reserve(ansatzPairs.size());
        std::transform(ansatzPairs.begin(),
                       ansatzPairs.end(),
                       std::back_inserter(connectivities),
                       [](const auto& rPair){return std::pair<Size,Size>(rPair.first, *rPair.second.begin());});
        _connectivityMap.emplace(iDim, std::move(connectivities));
    }

    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
template <concepts::OutputIterator<std::pair<Size,Size>> TOutputIt>
void AnsatzMap<TValue>::getPairs(BoundaryID boundary,
                                 TOutputIt itOutput) const noexcept
{
    const auto it = _connectivityMap.find(boundary.getDimension());
    if (it != _connectivityMap.end()) {
        if (boundary.getDirection()) {
            std::transform(it.value().begin(),
                           it.value().end(),
                           itOutput,
                           [](auto pair) -> std::pair<Size,Size>
                            {std::swap(pair.first, pair.second); return pair;});
        } else {
            std::copy(it.value().begin(),
                      it.value().end(),
                      itOutput);
        }
    }
}


template <class TValue>
Size AnsatzMap<TValue>::getPairCount(BoundaryID boundary) const noexcept
{
    const auto it = _connectivityMap.find(boundary.getDimension());
    if (it != _connectivityMap.end()) {
        return it->second.size();
    } else {
        return 0;
    }
}


template <maths::Expression TAnsatzSpace>
AnsatzMap<typename TAnsatzSpace::Value> makeAnsatzMap(Ref<const TAnsatzSpace> rAnsatzSpace,
                                                      std::span<const typename TAnsatzSpace::Value> samples,
                                                      utils::Comparison<typename TAnsatzSpace::Value> comparison)
{
    return AnsatzMap<typename TAnsatzSpace::Value>(rAnsatzSpace, samples, comparison);
}


} // namespace cie::fem


#endif
