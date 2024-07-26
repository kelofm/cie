#ifndef CIE_FEM_ASSEMBLER_IMPL_HPP
#define CIE_FEM_ASSEMBLER_IMPL_HPP

// --- External Includes ---
#include <tsl/robin_set.h>

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <queue>


namespace cie::fem {


template <class TVertexData,
          class TEdgeData,
          concepts::FunctionWithSignature<std::size_t,Ref<const typename Graph<TVertexData,TEdgeData>::Vertex>> TDoFCounter,
          concepts::FunctionWithSignature<void,Ref<const typename Graph<TVertexData,TEdgeData>::Edge>,Assembler::DoFPairIterator> TDoFPairFunctor>
void Assembler::addGraph(Ref<const Graph<TVertexData,TEdgeData>> rGraph,
                         TDoFCounter&& rDoFCounter,
                         TDoFPairFunctor&& rDoFMatcher)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Early exit if the graph is empty
    if (rGraph.empty()) {
        return;
    }

    using Vertex = typename Graph<TVertexData,TEdgeData>::Vertex;
    using Edge = typename Graph<TVertexData,TEdgeData>::Edge;

    std::queue<Ptr<const Vertex>> visitQueue {{&rGraph.vertices().front()}};
    tsl::robin_set<typename Vertex::ID> visited;
    DoFPairVector dofPairs;

    while (!visitQueue.empty()) {
        // Strip the next vertex to visit
        Ref<const Vertex> rVertex = *visitQueue.front();
        visited.insert(rVertex.id());
        visitQueue.pop();

        for (const auto edgeID : rVertex.edges()) {
            Ref<const Edge> rEdge = rGraph.find(edgeID).value();

            Ptr<const Vertex> pSource, pTarget;
            if (rEdge.source() == rVertex.id()) {
                pSource = &rVertex;
                pTarget = &rGraph.find(rEdge.target()).value();
                if (visited.emplace(rEdge.target()).second) {
                    visitQueue.push(pTarget);
                }
            } else {
                pSource = &rGraph.find(rEdge.source()).value();
                pTarget = &rVertex;
                if (visited.emplace(rEdge.source()).second) {
                    visitQueue.push(pSource);
                }
            }

            Ref<const Vertex> rSource = *pSource;
            Ref<const Vertex> rTarget = *pTarget;

            // Get DoF containers
            Ref<DoFMap::mapped_type> rSourceDoFs = _dofMap.emplace(rEdge.source(), DoFMap::mapped_type {}).first->second;
            Ref<DoFMap::mapped_type> rTargetDoFs = _dofMap.emplace(rEdge.target(), DoFMap::mapped_type {}).first->second;
            rSourceDoFs.resize(rDoFCounter(rSource));
            rTargetDoFs.resize(rDoFCounter(rTarget));

            // Get DoF connectivities
            dofPairs.clear();
            rDoFMatcher(rEdge, std::back_inserter(dofPairs));

            // Assign DoFs
            for (auto dofPair : dofPairs) {
                auto& rMaybeSourceDoF = rSourceDoFs[dofPair.first];
                auto& rMaybeTargetDoF = rTargetDoFs[dofPair.second];

                if (rMaybeSourceDoF.has_value()) {
                    if (rMaybeTargetDoF.has_value()) {
                        CIE_CHECK(
                            *rMaybeSourceDoF == *rMaybeTargetDoF,
                            "DoF assignment failure at edge " << rEdge.id()
                            << " between vertex " << rEdge.source() << " (local DoF " << dofPair.first << " assigned to global DoF " << *rMaybeSourceDoF << ")"
                            << " and vertex " << rEdge.target() << " (local DoF " << dofPair.second << " assigned to global DoF " << *rMaybeTargetDoF << ")"
                        );
                    } else {
                        rMaybeTargetDoF = rMaybeSourceDoF;
                    }
                } else if (rMaybeTargetDoF.has_value()) {
                    rMaybeSourceDoF = rMaybeTargetDoF;
                } else {
                    const auto iDoF = _dofCounter++;
                    rMaybeSourceDoF = iDoF;
                    rMaybeTargetDoF = iDoF;
                }
            } // for dofPair in dofPairs
        } // for rEdge in rVertex.edges()
    } // while visitQueue

    for (auto& rPair : _dofMap) for (auto& riDoF : rPair.second) if (!riDoF.has_value()) riDoF = _dofCounter++;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem


#endif
