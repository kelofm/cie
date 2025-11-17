#ifndef CIE_FEM_ASSEMBLER_IMPL_HPP
#define CIE_FEM_ASSEMBLER_IMPL_HPP

// --- External Includes ---
#include <tsl/robin_set.h>

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/Mutex.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <queue>
#include <algorithm> // lower_bound, sort, unique
#include <numeric> // inclusive_scan
#include <span>


namespace cie::fem {


template <class TVertexData,
          class TEdgeData,
          class TGraphData,
          concepts::FunctionWithSignature<std::size_t,Ref<const typename Graph<TVertexData,TEdgeData,TGraphData>::Vertex>> TDoFCounter,
          concepts::FunctionWithSignature<void,Ref<const typename Graph<TVertexData,TEdgeData,TGraphData>::Edge>,Assembler::DoFPairIterator> TDoFPairFunctor>
void Assembler::addGraph(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph,
                         TDoFCounter&& rDoFCounter,
                         TDoFPairFunctor&& rDoFMatcher)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Early exit if the graph is empty
    if (rGraph.empty()) {
        return;
    }

    using Vertex = typename Graph<TVertexData,TEdgeData,TGraphData>::Vertex;
    using Edge = typename Graph<TVertexData,TEdgeData,TGraphData>::Edge;

    std::queue<Ptr<const Vertex>> visitQueue {{&rGraph.vertices().front()}};
    tsl::robin_set<typename Vertex::ID> visited;
    DoFPairVector dofPairs;

    // Questionable reserve here.
    // On the one hand, it's extremely wasteful because it assumes that no
    // DoFs are shared between cells. On the other hand, it doesn't ensure
    // a final size because it assumes all vertices have an identical number
    // of DoFs.
    _dofMap.reserve(rGraph.vertices().size() * rDoFCounter(*visitQueue.front()));

    while (!visitQueue.empty()) {
        // Strip the next vertex to visit
        Ref<const Vertex> rVertex = *visitQueue.front();
        visited.insert(rVertex.id());
        visitQueue.pop();

        // Insert the vertex to handle disconnected cells.
        //_dofMap.emplace(rVertex.id(), DoFMap::mapped_type {}).first->second.resize(rDoFCounter(rVertex));
        _dofMap.emplace(rVertex.id(), DoFMap::mapped_type {}).first.value().resize(rDoFCounter(rVertex));

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
                rDoFs.resize(rDoFCounter(rSource));
                sourceDoFs = std::span<DoFMap::mapped_type::value_type>(rDoFs.data(), rDoFs.data() + rDoFs.size());
            }

            {
                Ref<DoFMap::mapped_type> rDoFs = _dofMap.emplace(rEdge.target(), DoFMap::mapped_type {}).first.value();
                rDoFs.resize(rDoFCounter(rTarget));
                targetDoFs = std::span<DoFMap::mapped_type::value_type>(rDoFs.data(), rDoFs.data() + rDoFs.size());
            }

            // Get DoF connectivities
            dofPairs.clear();
            rDoFMatcher(rEdge, std::back_inserter(dofPairs));

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

    //for (auto& rPair : _dofMap)
    //    for (auto& riDoF : rPair.second)
    //        if (!riDoF.has_value())
    //            riDoF = _dofCounter++;

    for (auto it=_dofMap.begin(); it!=_dofMap.end(); ++it)
        for (auto& riDoF : it.value())
            if (!riDoF.has_value())
                riDoF = _dofCounter++;

    CIE_END_EXCEPTION_TRACING
}


template <class TIndex, class TValue>
void Assembler::makeCSRMatrix(Ref<TIndex> rRowCount,
                              Ref<TIndex> rColumnCount,
                              Ref<DynamicArray<TIndex>> rRowExtents,
                              Ref<DynamicArray<TIndex>> rColumnIndices,
                              Ref<DynamicArray<TValue>> rNonzeros,
                              OptionalRef<mp::ThreadPoolBase> rThreadPool) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Initialize matrix containers
    rRowExtents.clear();
    rColumnIndices.clear();
    rNonzeros.clear();

    const TIndex rowCount = this->dofCount();
    rRowCount = rowCount;
    rColumnCount = rowCount;

    // Initialize mutexes
    // Each mutex is assigned a range of row indices it provides access to
    const TIndex threadCount = rThreadPool.has_value() ? rThreadPool.value().size() : 1;
    const TIndex minRowsPerMutex = rowCount / threadCount;
    const TIndex leftoverRowCount = rowCount % threadCount;
    DynamicArray<TIndex> mutexExtents;
    DynamicArray<mp::Mutex<tags::SMP>> mutexes(threadCount);

    mutexExtents.reserve(threadCount + 1);
    for (TIndex iThread=0; iThread<threadCount; ++iThread) {
        mutexExtents.push_back(mutexExtents.empty() ? 0 : mutexExtents.back()
                               + minRowsPerMutex
                               + (iThread < leftoverRowCount ? 1 : 0));
    }
    mutexExtents.push_back(rowCount);

    // Insert column indices into temporary row container
    DynamicArray<DynamicArray<TIndex>> columnIndices(rowCount);

    {
        const auto job = [&columnIndices, &mutexExtents, &mutexes](const auto& rIndices) -> void {
            for (const auto iRow : rIndices) {
                // Find which thread the row belongs to.
                const auto itMutexExtent = std::max(
                    mutexExtents.begin(),
                    std::lower_bound(mutexExtents.begin(),
                                    mutexExtents.end(),
                                    iRow) - 1
                );
                CIE_OUT_OF_RANGE_CHECK(itMutexExtent < mutexExtents.end());
                auto& rMutex = mutexes[std::distance(mutexExtents.begin(), itMutexExtent)];

                std::scoped_lock<mp::Mutex<tags::SMP>> lock(rMutex);
                CIE_OUT_OF_RANGE_CHECK(iRow < columnIndices.size());
                columnIndices[iRow].insert(columnIndices[iRow].end(), rIndices.begin(), rIndices.end());
            } // for iRow in rIndices
        };

        if (threadCount < 2) {
            for (const auto& rDofIndices : this->values()) job(rDofIndices);
        } else {
            const auto& rDofIndexContainers = this->values();
            mp::ParallelFor<>(rThreadPool.value())(rDofIndexContainers.begin(),
                                                   rDofIndexContainers.end(),
                                                   job);
        }
    }

    // Make column indices sorted and unique
    {
        const auto job = [] (Ref<DynamicArray<TIndex>> rIndices) {
            std::sort(rIndices.begin(), rIndices.end());
            rIndices.erase(std::unique(rIndices.begin(), rIndices.end()), rIndices.end());
        };

        if (threadCount < 2) {
            for (auto& rIndices : columnIndices) job(rIndices);
        } else {
            mp::ParallelFor<>(rThreadPool.value())(columnIndices, job);
        }
    }

    // Compute row extents
    rRowExtents.resize(rowCount + 1);
    rRowExtents.front() = 0;

    std::inclusive_scan(columnIndices.begin(),
                        columnIndices.end(),
                        rRowExtents.begin() + 1,
                        [](TIndex left, const auto& rRight) -> TIndex {return left + rRight.size();},
                        static_cast<TIndex>(0));

    // Copy column indices to CSR container
    rColumnIndices.resize(rRowExtents.back());
    rNonzeros.resize(rRowExtents.back());

    {
        const auto job = [&columnIndices, &rColumnIndices, &rRowExtents] (const TIndex iRow) {
            std::copy(columnIndices[iRow].begin(),
                      columnIndices[iRow].end(),
                      rColumnIndices.begin() + rRowExtents[iRow]);
        };

        if (threadCount < 2) {
            for (TIndex iRow=0; iRow<rowCount; ++iRow) job(iRow);
        } else {
            mp::ParallelFor<>(rThreadPool.value())(rowCount, job);
        }
    }

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem


#endif
