#pragma once

// --- External Includes ---
#include <tsl/robin_set.h>

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <queue>
#include <algorithm> // lower_bound, sort, unique
#include <span>


namespace cie::fem {


template <
    class TVertexData,
    class TEdgeData,
    class TGraphData,
    unsigned Dimension>
requires (CellLike<TVertexData> && CellBoundaryLike<TEdgeData>)
void Assembler::addGraph(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph,
                         Ref<const AnsatzMap<Dimension>> rAnsatzMap,
                         std::size_t dofsPerCell) {
    CIE_BEGIN_EXCEPTION_TRACING

    // Early exit if the graph is empty
    if (rGraph.empty()) {
        return;
    }

    using Vertex = typename Graph<TVertexData,TEdgeData,TGraphData>::Vertex;
    using Edge = typename Graph<TVertexData,TEdgeData,TGraphData>::Edge;

    std::queue<Ptr<const Vertex>> visitQueue {{&rGraph.vertices().front()}};
    tsl::robin_set<typename Vertex::ID> visited;
    DynamicArray<DofPair> dofPairs;

    // Questionable reserve here.
    // On the one hand, it's extremely wasteful because it assumes that no
    // DoFs are shared between cells. On the other hand, it doesn't ensure
    // a final size because it assumes all vertices have an identical number
    // of DoFs.
    _dofMap.reserve(rGraph.vertices().size() * dofsPerCell);

    while (!visitQueue.empty()) {
        // Strip the next vertex to visit
        Ref<const Vertex> rVertex = *visitQueue.front();
        visited.insert(rVertex.id());
        visitQueue.pop();

        // Insert the vertex to handle disconnected cells.
        _dofMap.emplace(rVertex.id(), DoFMap::mapped_type {}).first.value().resize(dofsPerCell);

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
            //Ref<DoFMap::mapped_type> rSourceDoFs = _dofMap.emplace(rEdge.source(), DoFMap::mapped_type {}).first->second;
            //Ref<DoFMap::mapped_type> rTargetDoFs = _dofMap.emplace(rEdge.target(), DoFMap::mapped_type {}).first->second;

            // Memory stability time.
            // OK, hear me out. I need to (potentially) insert two new items into the hash table (_dofMap),
            // which is not an iterator-stable operation. That said, the values this hash table stores are
            // dynamic arrays (i.e.: std::vector) that do not have small-array optimizations, which means
            // that the addresses of the objects they store remain stable unless a resizing is performed.
            // Thankfully, moving a dynamic array does not trigger a resizing so I can count on the contents
            // of my arrays to remain in place after inserting stuff into the hash table.
            // The catch is that even though most stdlib implementations define std::vector::iterator as
            // raw pointers (except for bool ...), it is not a requirement by the standard, nor is the
            // stability of iterators after a move. So raw pointers are the only way to go.
            std::span<DoFMap::mapped_type::value_type> sourceDoFs, targetDoFs;

            {
                Ref<DoFMap::mapped_type> rDoFs = _dofMap.emplace(
                    rEdge.source(),
                    DoFMap::mapped_type {}).first.value();
                rDoFs.resize(dofsPerCell);
                sourceDoFs = std::span<DoFMap::mapped_type::value_type>(rDoFs.data(), rDoFs.data() + rDoFs.size());
            }

            {
                Ref<DoFMap::mapped_type> rDoFs = _dofMap.emplace(rEdge.target(), DoFMap::mapped_type {}).first.value();
                rDoFs.resize(dofsPerCell);
                targetDoFs = std::span<DoFMap::mapped_type::value_type>(rDoFs.data(), rDoFs.data() + rDoFs.size());
            }

            // Get DoF connectivities.
            const OrientedBoundary<Dimension> sourceBoundary(rSource.data().axes(), rEdge.data().boundary());
            const OrientedBoundary<Dimension> targetBoundary(rTarget.data().axes(), rEdge.data().boundary());

            const auto itDofPairs = rAnsatzMap.findPairs(sourceBoundary, targetBoundary);
            dofPairs.resize(rAnsatzMap.pairCount(itDofPairs));
            rAnsatzMap.getPairs(itDofPairs, dofPairs);

            // Assign DoFs
            for (auto dofPair : dofPairs) {
                auto& rMaybeSourceDoF = sourceDoFs[dofPair.first];
                auto& rMaybeTargetDoF = targetDoFs[dofPair.second];

                if (rMaybeSourceDoF.has_value()) {
                    if (rMaybeTargetDoF.has_value()) {
                        CIE_CHECK(
                            *rMaybeSourceDoF == *rMaybeTargetDoF,
                            "DoF assignment failure at edge " << rEdge.id()
                            << " between vertex " << rEdge.source() << " (local DoF " << dofPair.first << " assigned to global DoF " << *rMaybeSourceDoF << ")"
                            << " and vertex " << rEdge.target() << " (local DoF " << dofPair.second << " assigned to global DoF " << *rMaybeTargetDoF << ")."
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

    for (auto it=_dofMap.begin(); it!=_dofMap.end(); ++it)
        for (auto& riDoF : it.value())
            if (!riDoF.has_value())
                riDoF = _dofCounter++;

    CIE_END_EXCEPTION_TRACING
}


template <
    TagLike TParallelism,
    concepts::Integer TIndex,
    concepts::Numeric TLocalScalar,
    concepts::Numeric TGlobalScalar>
void Assembler::addContribution(
    std::span<const TLocalScalar> contribution,
    VertexID cellID,
    std::span<const TIndex> rowExtents,
    std::span<const TIndex> columnIndices,
    std::span<TGlobalScalar> entries) const {
        static_assert(TParallelism::id() == tags::Serial::id() || TParallelism::id() == tags::SMP::id());
        const auto& rDofMap = this->operator[](cellID);
        const unsigned localSystemSize = rDofMap.size();
        for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
            for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
                const auto iRowBegin = rowExtents[rDofMap[iLocalRow]];
                const auto iRowEnd = rowExtents[rDofMap[iLocalRow] + 1];
                const auto itColumnIndex = std::lower_bound(
                    columnIndices.begin() + iRowBegin,
                    columnIndices.begin() + iRowEnd,
                    rDofMap[iLocalColumn]);
                CIE_OUT_OF_RANGE_CHECK(
                    itColumnIndex != columnIndices.begin() + iRowEnd
                    && *itColumnIndex == static_cast<TIndex>(rDofMap[iLocalColumn]));
                const auto iEntry = std::distance(columnIndices.begin(), itColumnIndex);

                if constexpr (TParallelism::id() == tags::SMP::id()) {
                    std::atomic_ref<TGlobalScalar>(entries[iEntry]) += contribution[iLocalRow * localSystemSize + iLocalColumn];
                } else if constexpr (TParallelism::id() == tags::Serial::id()) {
                    entries[iEntry] += contribution[iLocalRow * localSystemSize + iLocalColumn];
                } else {
                    static_assert(
                        std::is_same_v<TParallelism,void>,
                        "unsupported parallelism");
                }
            } // for iLocalColumn in range(ansatzBuffer.size)
        } // for iLocalRow in range(ansatzBuffer.size)
}


template <
    TagLike TParallelism,
    concepts::Numeric TLocalScalar,
    concepts::Numeric TGlobalScalar>
void Assembler::addContribution(
    std::span<const TLocalScalar> contribution,
    VertexID cellID,
    std::span<TGlobalScalar> entries) const
{
    const auto& rDofMap = this->operator[](cellID);
    for (unsigned iComponent=0u; iComponent<contribution.size(); ++iComponent) {
        const auto iRow = rDofMap[iComponent];

        if constexpr (TParallelism::id() == tags::SMP::id()) {
            std::atomic_ref<TGlobalScalar>(entries[iRow]) += contribution[iComponent];
        } else if constexpr (TParallelism::id() == tags::Serial::id()) {
            entries[iRow] += contribution[iComponent];
        } else {
            static_assert(
                std::is_same_v<TParallelism,void>,
                "unsupported parallelism");
        }
    }
}


} // namespace cie::fem
