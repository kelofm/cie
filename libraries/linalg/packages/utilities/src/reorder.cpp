// --- External Includes ---
#include <tsl/robin_set.h>

// --- Linalg Includes ---
#include "packages/utilities/inc/reorder.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <format>
#include <algorithm>
#include <vector>
#include <deque>


namespace cie {


template <concepts::Integer TIndex, concepts::Numeric TValue>
void makeReordering(
    std::span<TIndex> map,
    std::span<const TIndex> rowExtents,
    std::span<const TIndex> columnIndices,
    std::span<const TValue> entries,
    ReorderingStrategy strategy,
    OptionalRef<mp::ThreadPoolBase> rMaybePool) {
        // Sanity checks.
        CIE_CHECK(
            !rowExtents.empty(),
            "row extents of a CSR matrix must have at least 1 entry")
        CIE_CHECK(
            map.size() == rowExtents.size() - 1,
            std::format(
                "reordering map size ({}) does not match the number of rows ({}) in the provided CSR matrix",
                map.size(),
                rowExtents.size() - 1))
        CIE_CHECK(
            std::is_sorted(rowExtents.begin(), rowExtents.end()),
            "row extents of a CSR matrix must be sorted")
        CIE_CHECK(
            columnIndices.size() == entries.size(),
            std::format(
                "mismatched number of column indices ({}) and entries ({}) in CSR matrix",
                columnIndices.size(),
                entries.size()))

        if (columnIndices.empty()) return;
        CIE_CHECK(
            static_cast<std::size_t>(rowExtents.back()) == columnIndices.size(),
            std::format(
                "the last row's sentinel ({}) of the CSR matrix do not match the number of entries ({})",
                rowExtents.back(),
                columnIndices.size()))

        if (strategy == ReorderingStrategy::None) return;

        CIE_BEGIN_EXCEPTION_TRACING

        const std::size_t rowCount = rowExtents.size() - 1ul;
        tsl::robin_set<TIndex> visitedOrQueued;
        std::vector<TIndex> vertexDegrees(rowCount);
        std::deque<TIndex> toVisit;

        auto itMap = map.begin();

        // Compute vertex degrees and find the minimum one.
        {
            TIndex minDegree = std::numeric_limits<TIndex>::max();
            TIndex iMinRow = 0ul;

            if (rMaybePool.has_value()) {
                auto loop = mp::ParallelFor<TIndex>(rMaybePool.value())
                    .firstPrivate(minDegree, iMinRow);
                loop(
                    rowCount,
                    [rowExtents, &vertexDegrees] (TIndex iRow, Ref<TIndex> rMinDegree, Ref<TIndex> riMinRow) {
                        const TIndex rowSize = rowExtents[iRow + 1] - rowExtents[iRow];
                        vertexDegrees[iRow] = rowSize;
                        if (rowSize < rMinDegree) {
                            rMinDegree = rowSize;
                            riMinRow = iRow;
                        }
                    });

                for (const auto& rThreadStorage : loop.getPool().getStorage()) {
                    if (std::get<0>(rThreadStorage.values()) < minDegree) {
                        minDegree = std::get<0>(rThreadStorage.values());
                        iMinRow = std::get<1>(rThreadStorage.values());
                    }
                } // for rThreadStorage
            } /*if rMaybePool*/ else {
                for (TIndex iRow=0; iRow<static_cast<TIndex>(rowCount); ++iRow) {
                    const TIndex rowSize = rowExtents[iRow + 1] - rowExtents[iRow];
                    vertexDegrees[iRow] = rowSize;
                    if (rowSize < minDegree) {
                        minDegree = rowSize;
                        iMinRow = iRow;
                    }
                } // for iRow in range(rowCount)
            }

            // Begin traversal at the vertex with minimum degree.
            toVisit.push_back(iMinRow);
            visitedOrQueued.insert(iMinRow);
        }

        while (!toVisit.empty() || visitedOrQueued.size() < rowCount) {
            // Strip the next vertex to visit.
            const TIndex iRow = toVisit.front();
            toVisit.pop_front();

            // Insert unvisited neighbors into the visit queue
            // in ascending order of their vertex degree.
            const TIndex iEntryBegin = rowExtents[iRow];
            const TIndex iEntryEnd = rowExtents[iRow + 1];
            const std::size_t sizeBeforeInsertion = toVisit.size();

            for (TIndex iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                const TIndex iVertex = columnIndices[iEntry];
                if (visitedOrQueued.find(iVertex) == visitedOrQueued.end()) {
                    toVisit.push_back(iVertex);
                    visitedOrQueued.insert(iVertex);
                }
            }

            std::sort(
                toVisit.begin() + sizeBeforeInsertion,
                toVisit.end(),
                [&vertexDegrees](TIndex left, TIndex right) -> bool {
                    return vertexDegrees[left] < vertexDegrees[right];
                });

            *itMap++ = iRow;

            // If the visit queue is empty, find the unvisited vertex
            // with the lowest degree and add it to the queue.
            if (toVisit.empty()) {
                std::pair<TIndex,TIndex> minElement {
                    std::numeric_limits<TIndex>::max(),
                    std::numeric_limits<TIndex>::max()};

                for (TIndex iRow=0; static_cast<std::size_t>(iRow)<rowCount; ++iRow) {
                    if (vertexDegrees[iRow] < minElement.first) {
                        if (visitedOrQueued.find(iRow) == visitedOrQueued.end()) {
                            minElement = std::make_pair(vertexDegrees[iRow], iRow);
                        }
                    }
                } // for iRow in range(rowCount)

                if (minElement.second != std::numeric_limits<TIndex>::max()) {
                    toVisit.push_back(minElement.second);
                    visitedOrQueued.insert(minElement.second);
                }
            } // if toVisit.empty()
        } // while visitedOrQueued.size() < rowCount || !toVisit.empty()

        CIE_CHECK(
            itMap == map.end(),
            std::format(
                "failed to build index map during reordering (reordererd {}/{} rows)",
                std::distance(map.begin(), itMap),
                rowCount))

        if (strategy == ReorderingStrategy::ReverseCuthillMcKee) {
            std::reverse(map.begin(), map.end());
        }

        CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, concepts::Numeric TValue, bool Reverse>
void reorder(
    std::span<const TIndex> map,
    std::span<TIndex> rowExtents,
    std::span<TIndex> columnIndices,
    std::span<TValue> entries,
    [[maybe_unused]] OptionalRef<mp::ThreadPoolBase> rMaybePool) {
        // Sanity checks.
        CIE_CHECK(
            !rowExtents.empty(),
            "row extents of a CSR matrix must have at least 1 entry")
        CIE_CHECK(
            std::is_sorted(rowExtents.begin(), rowExtents.end()),
            "row extents of a CSR matrix must be sorted")
        CIE_CHECK(
            columnIndices.size() == entries.size(),
            std::format(
                "mismatched number of column indices ({}) and entries ({}) in CSR matrix",
                columnIndices.size(),
                entries.size()))
        CIE_CHECK(
            map.size() == rowExtents.size() - 1,
            std::format(
                "map size ({}) does not match the number of rows ({}) in the CSR matrix",
                map.size(),
                rowExtents.size() - 1))

        if (columnIndices.empty()) return;
        CIE_CHECK(
            static_cast<std::size_t>(rowExtents.back()) == columnIndices.size(),
            std::format(
                "the last row's sentinel ({}) of the CSR matrix do not match the number of entries ({})",
                rowExtents.back(),
                columnIndices.size()))

        CIE_BEGIN_EXCEPTION_TRACING

        const std::size_t rowCount = rowExtents.size() - 1;

        std::vector<TIndex> reverseMap, newRowExtents(rowExtents.size()), newColumnIndices(columnIndices.size());
        std::vector<TValue> newEntries(entries.size());
        reverseMap.reserve(map.size());
        reverseMap.front() = static_cast<TIndex>(0);

        // Compute the reverse map.
        const auto reverseMapKernel = [&map, &reverseMap] (TIndex iRow) -> void {
            const TIndex iOldRow = map[iRow];
            reverseMap[iOldRow] = iRow;
        };

        if (rMaybePool.has_value()) {
            mp::ParallelFor<TIndex>(rMaybePool.value()).operator()(
                rowCount,
                reverseMapKernel);
        } /*if rMaybePool*/ else {
            for (TIndex iRow=0; iRow<static_cast<TIndex>(rowCount); ++iRow) reverseMapKernel(iRow);
        }

        const std::span<const TIndex> mapView = Reverse ? reverseMap : map;
        const std::span<const TIndex> reverseMapView = Reverse ? map : reverseMap;

        // Compute new row extents.
        for (TIndex iRow=0; iRow<static_cast<TIndex>(rowCount); ++iRow) {
            const TIndex iOldRow = mapView[iRow];
            const TIndex rowSize = rowExtents[iOldRow + 1] - rowExtents[iOldRow];
            newRowExtents[iRow + 1] = newRowExtents[iRow] + rowSize;
        } // for iRow in range(rowCount)

        // Define the reordering kernel.
        const auto reorderingKernel = [&] (TIndex iRow) -> void {
            const TIndex iOldRow = mapView[iRow];

            const TIndex iOldEntryBegin = rowExtents[iOldRow];
            const TIndex iNewEntryBegin = newRowExtents[iRow];
            const TIndex iOldEntryEnd = rowExtents[iOldRow + 1];

            for (TIndex iOldEntry=iOldEntryBegin, iNewEntry=iNewEntryBegin; iOldEntry<iOldEntryEnd; ++iOldEntry, ++iNewEntry) {
                newColumnIndices[iNewEntry] = reverseMapView[columnIndices[iOldEntry]];
                newEntries[iNewEntry] = entries[iOldEntry];
            } // for iOldEntry, iNewEntry
        }; // reorderingKernel

        // Execute the reordering kernel.
        if (rMaybePool.has_value()) {
            mp::ParallelFor<TIndex>(rMaybePool.value()).operator()(
                rowCount,
                reorderingKernel);
        } /*if rMaybePool.has_value()*/ else {
            for (TIndex iRow=0; iRow<static_cast<TIndex>(rowCount); ++iRow) reorderingKernel(iRow);
        }

        // Overwrite input arrays.
        std::copy(
            newRowExtents.begin(),
            newRowExtents.end(),
            rowExtents.begin());
        std::copy(
            newColumnIndices.begin(),
            newColumnIndices.end(),
            columnIndices.begin());
        std::copy(
            newEntries.begin(),
            newEntries.end(),
            entries.begin());

        CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, concepts::Numeric TValue, bool Reverse>
void reorder(
    std::span<const TIndex> map,
    std::span<TValue> array,
    OptionalRef<mp::ThreadPoolBase> rMaybePool) {
        CIE_BEGIN_EXCEPTION_TRACING
        std::vector<TValue> output(array.size());
        const auto kernel = [map, array, &output] (TIndex iEntry) -> void {
            if constexpr (Reverse) {
                output[map[iEntry]] = array[iEntry];
            } else {
                output[iEntry] = array[map[iEntry]];
            }
        };

        if (rMaybePool.has_value()) {
            mp::ParallelFor<TIndex>(rMaybePool.value()).operator()(array.size(), kernel);
        } else {
            for (TIndex iEntry=0; iEntry<static_cast<TIndex>(array.size()); ++iEntry) kernel(iEntry);
        }

        std::copy(
            output.begin(),
            output.end(),
            array.begin());
        CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, concepts::Numeric TValue>
void reverseReorder(
    std::span<const TIndex> map,
    std::span<TIndex> rowExtents,
    std::span<TIndex> columnIndices,
    std::span<TValue> entries,
    OptionalRef<mp::ThreadPoolBase> rMaybePool) {
        CIE_BEGIN_EXCEPTION_TRACING
        reorder<TIndex,TValue,true>(
            map,
            rowExtents,
            columnIndices,
            entries,
            rMaybePool);
        CIE_END_EXCEPTION_TRACING
}


template <concepts::Integer TIndex, concepts::Numeric TValue>
void reverseReorder(
    std::span<const TIndex> map,
    std::span<TValue> array,
    OptionalRef<mp::ThreadPoolBase> rMaybePool) {
        CIE_BEGIN_EXCEPTION_TRACING
        reorder<TIndex,TValue,true>(
            map,
            array,
            rMaybePool);
        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_REORDER(TIndex, TValue)         \
    template void makeReordering<TIndex,TValue>(        \
        std::span<TIndex>,                              \
        std::span<const TIndex>,                        \
        std::span<const TIndex>,                        \
        std::span<const TValue>,                        \
        ReorderingStrategy,                             \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reorder<TIndex,TValue,false>(         \
        std::span<const TIndex>,                        \
        std::span<TIndex>,                              \
        std::span<TIndex>,                              \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reorder<TIndex,TValue,true>(          \
        std::span<const TIndex>,                        \
        std::span<TIndex>,                              \
        std::span<TIndex>,                              \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reverseReorder<TIndex,TValue>(        \
        std::span<const TIndex>,                        \
        std::span<TIndex>,                              \
        std::span<TIndex>,                              \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reorder<TIndex,TValue,false>(         \
        std::span<const TIndex>,                        \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reorder<TIndex,TValue,true>(          \
        std::span<const TIndex>,                        \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);               \
    template void reverseReorder<TIndex,TValue>(        \
        std::span<const TIndex>,                        \
        std::span<TValue>,                              \
        OptionalRef<mp::ThreadPoolBase>);

CIE_INSTANTIATE_REORDER(int, float)
CIE_INSTANTIATE_REORDER(int, double)
CIE_INSTANTIATE_REORDER(std::size_t, float)
CIE_INSTANTIATE_REORDER(std::size_t, double)

#undef CIE_INSTANTIATE_REORDER

} // namespace cie
