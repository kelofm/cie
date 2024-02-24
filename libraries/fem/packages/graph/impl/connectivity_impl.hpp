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
void scanConnectivities(Ref<const TAnsatzSpace> r_ansatzSpace,
                        TFunctor&& r_functor,
                        Ptr<const typename TAnsatzSpace::Value> p_sampleBegin,
                        Ptr<const typename TAnsatzSpace::Value> p_sampleEnd,
                        typename TAnsatzSpace::Value tolerance)
{
    CIE_BEGIN_EXCEPTION_TRACING
    using Value = typename TAnsatzSpace::Value;
    constexpr unsigned Dimension = TAnsatzSpace::Dimension;
    static_assert(0 < Dimension);

    const Size ansatzSize = r_ansatzSpace.size();
    const Size numberOfSamples = std::distance(p_sampleBegin, p_sampleEnd);

    StaticArray<Value,Dimension> argument;
    StaticArray<Size,Dimension-1> argumentState;
    DynamicArray<Value> valueBuffer(ansatzSize);

    std::fill(argumentState.begin(),
              argumentState.end(),
              static_cast<Size>(0));

    constexpr unsigned maxBoundaries = 2 * Dimension;
    BoundaryID boundaryID;
    tsl::robin_set<std::pair<BoundaryID,Size>> activeBoundaries;

    for (unsigned i_boundary=0; i_boundary<maxBoundaries; ++i_boundary, ++boundaryID) {
        activeBoundaries.clear();
        const unsigned i_dim = boundaryID.getDimension();

        do {
            // Compute the current sample point
            Size i_state = 0;
            for (unsigned i=0; i<Dimension; ++i) {
                if (i == i_dim) {
                    argument[i] = static_cast<Value>(boundaryID.getDirection() ? 1 : -1);
                } else {
                    argument[i] = *(p_sampleBegin + argumentState[i_state++]);
                }
            } // for i in range(Dimension)

            // Evaluate the ansatz space
            r_ansatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data());

            for (Size i_ansatz=0; i_ansatz<ansatzSize; ++i_ansatz) {
                if (tolerance < valueBuffer[i_ansatz]) {
                    activeBoundaries.emplace(boundaryID, i_ansatz);
                }
            }
        } while (maths::CartesianProduct<Dimension-1>::next(numberOfSamples, argumentState.data()));

        for (auto pair : activeBoundaries) {
            r_functor(pair.first, pair.second);
        }
    } // for i_boundary in range(maxBoundaries)
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
template <maths::Expression TAnsatzSpace>
requires (std::is_same_v<typename TAnsatzSpace::Value,TValue>)
AnsatzMap<TValue>::AnsatzMap(Ref<const TAnsatzSpace> r_ansatzSpace,
                             std::span<const TValue> samples,
                             utils::Comparison<TValue> comparison)
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr unsigned Dimension = TAnsatzSpace::Dimension;
    static_assert(0 < Dimension);

    using Value = typename TAnsatzSpace::Value;
    const Size ansatzSize = r_ansatzSpace.size();   // <== total number of ansatz functions
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
    for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
        // Assume every ansatz function vanishes on both boundaries
        nonzeros.clear();

        // Assume every ansatz function is coincident with every other function
        // (consistent with the assumption of all functions vanishing on both boundaries)
        {
            tsl::robin_set<Size> all;
            all.reserve(ansatzSize);
            for (Size i=0; i<ansatzSize; ++i) all.insert(i);
            for (Size i_ansatz=0; i_ansatz<ansatzSize; ++i_ansatz) ansatzPairs.insert_or_assign(i_ansatz, all);
        }

        do {
            // Compute the current sample point
            unsigned i_state = 0;
            for (unsigned i=0; i<Dimension; ++i) {
                if (i != i_dim) {
                    argument[i] = samples[argumentState[i_state++]];
                }
            } // for i in range(Dimension)

            // Evaluate the ansatz space on the negative boundary
            argument[i_dim] = Value(-1);
            r_ansatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data());

            // Evaluate the ansatz space on the positive boundary
            argument[i_dim] = Value(1);
            r_ansatzSpace.evaluate(argument.begin(),
                                   argument.end(),
                                   valueBuffer.data() + ansatzSize);

            // Compare ansatz function values at opposite boundaries
            for (unsigned i_ansatz=0; i_ansatz<ansatzSize; ++i_ansatz) {
                // Check whether the negative boundary is nonzero
                if (!comparison.equal(valueBuffer[i_ansatz], Value(0))) {
                    nonzeros.emplace(BoundaryID(i_dim, false), i_ansatz);
                }

                // Check whether the positive boundary is nonzero
                if (!comparison.equal(valueBuffer[i_ansatz + ansatzSize], Value(0))) {
                    nonzeros.emplace(BoundaryID(i_dim, true), i_ansatz);
                }
            }

            // Update the coincidence map
            for (auto it=ansatzPairs.begin(); it!=ansatzPairs.end(); ++it) {
                const Size i_negative = it.key(); // <== index of the ansatz function on the negative boundary
                DynamicArray<Size> erase;

                for (Size i_positive : it.value()) {
                    if (!comparison.equal(valueBuffer[i_negative], valueBuffer[i_positive + ansatzSize])) {
                        erase.push_back(i_positive);
                    } // if not coincident
                } // for i_positive in pair.second

                for (auto i_erase : erase) {
                    it.value().erase(i_erase);
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
        for (Size i_ansatz=0; i_ansatz<ansatzSize; ++i_ansatz) {
            if (nonzeros.find({{i_dim, false}, i_ansatz}) == nonzeros.end()) {
                ansatzPairs.erase(i_ansatz);
            }
        }

        // Check whether each ansatz function that doesn't vanish on the negative
        // boundary has exactly ONE coincident ansatz function on the positive boundary.
        for (auto it=ansatzPairs.begin(); it!=ansatzPairs.end(); ++it) {
            if (it.value().size() != 1) {
                std::stringstream message;
                message << "Ansatz function " << it.key() << " at the negative boundary of dimension " << i_dim
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
        for (auto [boundary, i_ansatz] : nonzeros) {
            if (boundary.getDirection()) {
                bool foundCoincidentPair = false;
                for (const auto& r_pair : ansatzPairs) {
                    if (r_pair.second.find(i_ansatz) != r_pair.second.end()) {
                        foundCoincidentPair = true;
                        break;
                    }
                }
                if (!foundCoincidentPair) {
                    std::stringstream message;
                    message << "Ansatz function " << i_ansatz
                            << " does not vanish on the positive boundary of dimension " << i_dim
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
                       [](const auto& r_pair){return std::pair<Size,Size>(r_pair.first, *r_pair.second.begin());});
        _connectivityMap.emplace(i_dim, std::move(connectivities));
    }

    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
template <concepts::OutputIterator<std::pair<Size,Size>> TOutputIt>
void AnsatzMap<TValue>::getPairs(BoundaryID boundary,
                                 TOutputIt it_output) const noexcept
{
    const auto it = _connectivityMap.find(boundary.getDimension());
    if (it != _connectivityMap.end()) {
        if (boundary.getDirection()) {
            std::transform(it.value().begin(),
                           it.value().end(),
                           it_output,
                           [](auto pair) -> std::pair<Size,Size>
                            {std::swap(pair.first, pair.second); return pair;});
        } else {
            std::copy(it.value().begin(),
                      it.value().end(),
                      it_output);
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
AnsatzMap<typename TAnsatzSpace::Value> makeAnsatzMap(Ref<const TAnsatzSpace> r_ansatzSpace,
                                                      std::span<const typename TAnsatzSpace::Value> samples,
                                                      utils::Comparison<typename TAnsatzSpace::Value> comparison)
{
    return AnsatzMap<typename TAnsatzSpace::Value>(r_ansatzSpace, samples, comparison);
}


} // namespace cie::fem


#endif
