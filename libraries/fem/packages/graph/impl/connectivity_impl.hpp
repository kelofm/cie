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


namespace std {


inline bool operator==(std::pair<cie::fem::BoundaryID,unsigned> left,
                       std::pair<cie::fem::BoundaryID,unsigned> right) noexcept
{
    return (left.first == right.first) && (left.second == right.second);
}


inline bool operator!=(std::pair<cie::fem::BoundaryID,unsigned> left,
                       std::pair<cie::fem::BoundaryID,unsigned> right) noexcept
{
    return left.first != right.first && left.second != right.second;
}


template <>
struct hash<std::pair<cie::fem::BoundaryID,unsigned>>
{
    auto operator()(std::pair<cie::fem::BoundaryID,unsigned> item) const
    {
        const auto tmp = hash<cie::fem::BoundaryID>()(item.first);
        return tmp ^ (hash<unsigned>()(item.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
    }
};


} // namespace std


namespace cie::fem {


template <class TScalarExpression, unsigned Dimension, concepts::CallableWith<BoundaryID,unsigned> TFunctor>
void scanConnectivities(Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace,
                        TFunctor&& r_functor,
                        Ptr<const typename TScalarExpression::Value> p_sampleBegin,
                        Ptr<const typename TScalarExpression::Value> p_sampleEnd,
                        typename TScalarExpression::Value tolerance)
{
    static_assert(0 < Dimension);

    using Value = typename TScalarExpression::Value;
    const unsigned ansatzSize = r_ansatzSpace.size();
    const unsigned numberOfSamples = std::distance(p_sampleBegin, p_sampleEnd);

    StaticArray<Value,Dimension> argument;
    StaticArray<unsigned,Dimension-1> argumentState;
    DynamicArray<Value> valueBuffer(ansatzSize);

    std::fill(argumentState.begin(),
              argumentState.end(),
              static_cast<unsigned>(0));

    constexpr unsigned maxBoundaries = 2 * Dimension;
    BoundaryID boundaryID;
    tsl::robin_set<std::pair<BoundaryID,unsigned>> activeBoundaries;

    for (unsigned i_boundary=0; i_boundary<maxBoundaries; ++i_boundary, ++boundaryID) {
        activeBoundaries.clear();
        const unsigned i_dim = boundaryID.getDimension();

        do {
            // Compute the current sample point
            unsigned i_state = 0;
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

            for (unsigned i_ansatz=0; i_ansatz<ansatzSize; ++i_ansatz) {
                if (tolerance < valueBuffer[i_ansatz]) {
                    activeBoundaries.emplace(boundaryID, i_ansatz);
                }
            }
        } while (maths::CartesianProduct<Dimension-1>::next(numberOfSamples, argumentState.data()));

        for (auto pair : activeBoundaries) {
            r_functor(pair.first, pair.second);
        }
    } // for i_boundary in range(maxBoundaries)
}


} // namespace cie::fem


#endif
