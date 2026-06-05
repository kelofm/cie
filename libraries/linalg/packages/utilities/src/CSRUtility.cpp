// --- External Includes ---
#include "tsl/robin_set.h"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Linalg Includes ---
#include "packages/utilities/inc/CSRUtility.hpp"

// --- STL Includes ---
#include <format>
#include <algorithm>
#include <numeric>
#include <iostream>


namespace cie::linalg {


template <class T, class I>
CSRUtility<T,I>::CSRUtility(OptionalRef<mp::ThreadPoolBase> rMaybeThreads) noexcept
    : _rMaybeThreads(rMaybeThreads)
{}


template <class T, class I>
void CSRUtility<T,I>::sum(
    CSRView<T,const I> left,
    CSRView<const T,const I> right,
    T scale) {
        // Sanity checks.
        CIE_CHECK(
            left.rowCount() == right.rowCount()
            && left.columnCount() == right.columnCount(),
            std::format(
                "incompatible matrix sum [{}x{}] += [{}x{}]",
                left.rowCount(), left.columnCount(),
                right.rowCount(), right.columnCount()))

        CIE_CHECK(
            right.entries().size() <= left.entries().size(),
            "right operand to matrix sum has more entries than the left one")

        if (right.entries().empty()) return;

        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<const I> pLeftColumnBegin = left.columnIndices().data();
            Ptr<T> pLeftMutableEntryBegin  = left.entries().data();

            const auto job = [left, right, pLeftColumnBegin, pLeftMutableEntryBegin] (I iRow, const auto& op) -> void {
                const I iLeftEntryBegin        = left.rowExtents()[iRow];
                const I iLeftEntryEnd          = left.rowExtents()[iRow + 1];
                const Ptr<const I> itLeftColumnBegin = pLeftColumnBegin + iLeftEntryBegin;
                const Ptr<const I> itLeftColumnEnd   = pLeftColumnBegin + iLeftEntryEnd;

                const I iRightEntryBegin = right.rowExtents()[iRow];
                const I iRightEntryEnd   = right.rowExtents()[iRow + 1];

                Ptr<const I> itLeftColumn = itLeftColumnBegin;
                for (I iRightEntry=iRightEntryBegin; iRightEntry<iRightEntryEnd; ++iRightEntry) {
                    const I iColumn = right.columnIndices()[iRightEntry];
                    itLeftColumn = std::lower_bound(
                        itLeftColumn,
                        itLeftColumnEnd,
                        iColumn);
                    CIE_CHECK(
                        itLeftColumn != itLeftColumnEnd && (*itLeftColumn) == iColumn,
                        std::format(
                            "left hand operand in matrix sum has no entry in row {} at column {}",
                            iRow, iColumn))
                    const I iLeftEntry = static_cast<I>(std::distance<Ptr<const I>>(
                        left.columnIndices().data(),
                        itLeftColumn));
                    pLeftMutableEntryBegin[iLeftEntry] += op(right.entries()[iRightEntry]);
                } // for iRightEntry in range(iRightEntryBegin, iRightEntryEnd)
            }; // op

            if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                if (scale == static_cast<T>(1)) {
                    const auto op = [] (T r) -> T {return r;};
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        left.rowCount(),
                        [job, op] (I iRow) {job(iRow, op);});
                } else {
                    const auto op = [scale] (T r) -> T {return scale * r;};
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        left.rowCount(),
                        [job, op] (I iRow) {job(iRow, op);});
                }
            } else {
                if (scale == static_cast<T>(1)) {
                    const auto op = [] (T r) -> T {return r;};
                    for (I iRow=0; iRow<left.rowCount(); ++iRow)
                        job(iRow, op);
                } else {
                    const auto op = [scale] (T r) -> T {return scale * r;};
                    for (I iRow=0; iRow<left.rowCount(); ++iRow)
                        job(iRow, op);
                }

            }
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void CSRUtility<T,I>::sum(
    Ref<std::vector<I>> rTargetRowExtents,
    Ref<std::vector<I>> rTargetColumnIndices,
    Ref<std::vector<T>> rTargetEntries,
    CSRView<const T,const I> right,
    T targetScale,
    T rightScale) {
        CIE_BEGIN_EXCEPTION_TRACING
            CIE_CHECK(!rTargetRowExtents.empty(), "")
            const I rowCount = rTargetRowExtents.size() - 1;

            // Sanity checks.
            CIE_CHECK(
                rowCount == right.columnCount(),
                std::format(
                    "provided incompatible matrices to sum: ({}x{}) + ({}x{})",
                    rowCount, "?",
                    right.rowCount(), right.columnCount()))
            if (!right.columnCount() || !right.rowCount() || !right.entries().size()) return;

            std::vector<I> swapRowExtents(rTargetRowExtents.size(), I(0)), swapColumnIndices;
            std::vector<T> swapEntries;

            // Merge rows into separate containers.
            {
                std::vector<std::vector<std::pair<I,T>>> rows(rowCount);

                const auto sumJob = [&] (I iRow) -> void {
                    const I iLeftBegin  = rTargetRowExtents[iRow];
                    const I iLeftEnd    = rTargetRowExtents[iRow + 1];
                    const I iRightBegin = right.rowExtents()[iRow];
                    const I iRightEnd   = right.rowExtents()[iRow + 1];
                    auto& rRow = rows[iRow];
                    rRow.reserve((iLeftEnd - iLeftBegin) + (iRightEnd - iRightBegin));
                    for (I iEntry=iLeftBegin; iEntry<iLeftEnd; ++iEntry)
                        rRow.emplace_back(
                            rTargetColumnIndices[iEntry],
                            targetScale * rTargetEntries[iEntry]);
                    for (I iEntry=iRightBegin; iEntry<iRightEnd; ++iEntry)
                        rRow.emplace_back(
                            right.columnIndices()[iEntry],
                            rightScale * right.entries()[iEntry]);
                    std::sort(
                        rRow.begin(),
                        rRow.end(),
                        [] (const auto& rLeft, const auto& rRight) {return rLeft.first < rRight.first;});
                    auto itFirst = rRow.begin();
                    for (auto it=rRow.begin()+1; it<rRow.end(); ++it) {
                        if (itFirst->first == it->first) itFirst->second += it->second;
                        else (itFirst = it);
                    }
                    rRow.erase(
                        std::unique(
                            rRow.begin(),
                            rRow.end(),
                            [] (const auto& rLeft, const auto& rRight) {return rLeft.first == rRight.first;}),
                        rRow.end());
                }; // sumJob

                if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        rowCount,
                        sumJob);
                } else {
                    for (I iRow=0; iRow<rowCount; ++iRow)
                        sumJob(iRow);
                }

                // Compute new row extents.
                swapRowExtents.front() = I(0);
                for (I iRow=0; iRow<rowCount; ++iRow)
                    swapRowExtents[iRow + 1] = swapRowExtents[iRow] + rows[iRow].size();

                // Fill column indices and entries.
                swapColumnIndices.resize(swapRowExtents.back());
                swapEntries.resize(swapRowExtents.back());
                const auto copyJob = [&] (I iRow) -> void {
                    auto itColumnIndex = swapColumnIndices.data() + swapRowExtents[iRow];
                    auto itEntry       = swapEntries.data() + swapRowExtents[iRow];
                    for (auto itPair=rows[iRow].begin(); itPair!=rows[iRow].end(); ++itPair, ++itColumnIndex, ++itEntry) {
                        *itColumnIndex = itPair->first;
                        *itEntry       = itPair->second;
                    } // for itPair
                }; // copyJob

                if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        rowCount,
                        copyJob);
                } else {
                    for (I iRow=0; iRow<rowCount; ++iRow)
                        copyJob(iRow);
                }
            }

            // Construct the new matrix.
            rTargetRowExtents = std::move(swapRowExtents);
            rTargetColumnIndices = std::move(swapColumnIndices);
            rTargetEntries = std::move(swapEntries);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
bool CSRUtility<T,I>::hasSparseOverlap(
    std::span<const I> left,
    std::span<const I> right) {
        if (left.empty() || right.empty()) return false;
        Ptr<const I> itLeft = left.data();
        Ptr<const I> itRight = right.data();
        const Ptr<const I> itLeftEnd = left.data() + left.size();
        const Ptr<const I> itRightEnd = right.data() + right.size();
        for (; itLeft<itLeftEnd && itRight<itRightEnd;) {
            if (*itLeft == *itRight) return true;
            //else if (*itLeft < *itRight)
            //    itLeft = std::lower_bound(
            //        itLeft + 1,
            //        itLeftEnd,
            //        *itRight);
            //else
            //    itRight = std::lower_bound(
            //        itRight +  1,
            //        itRightEnd,
            //        *itLeft);
            else if (*itLeft < *itRight)
                ++itLeft;
            else
                ++itRight;
        }
        return false;
}


template <class T, class I>
T CSRUtility<T,I>::sparseInnerProduct(
    I iRow,
    I iColumn,
    Ref<const CSRView<const T, const I>> rMatrix) {
        T product = 0;
        const std::span<const I> leftIndices = {
            rMatrix.columnIndices().data() + rMatrix.rowExtents()[iRow],
            rMatrix.columnIndices().data() + rMatrix.rowExtents()[iRow + 1]};
        const std::span<const I> rightIndices = {
            rMatrix.columnIndices().data() + rMatrix.rowExtents()[iColumn],
            rMatrix.columnIndices().data() + rMatrix.rowExtents()[iColumn + 1]};
        if (leftIndices.empty() || rightIndices.empty()) return product;
        Ptr<const I> itLeft = leftIndices.data();
        Ptr<const I> itRight = rightIndices.data();
        const Ptr<const I> itLeftEnd = leftIndices.data() + leftIndices.size();
        const Ptr<const I> itRightEnd = rightIndices.data() + rightIndices.size();
        for (; itLeft<itLeftEnd && itRight<itRightEnd;) {
            if (*itLeft == *itRight) {
                const I iLeftEntry  = std::distance<Ptr<const I>>(rMatrix.columnIndices().data(), itLeft);
                const I iRightEntry = std::distance<Ptr<const I>>(rMatrix.columnIndices().data(), itRight);
                product += rMatrix.entries()[iLeftEntry] * rMatrix.entries()[iRightEntry];
                ++itLeft;
                ++itRight;
            }
            //else if (*itLeft < *itRight)
            //    itLeft = std::lower_bound(
            //        itLeft + 1,
            //        itLeftEnd,
            //        *itRight);
            //else
            //    itRight = std::lower_bound(
            //        itRight +  1,
            //        itRightEnd,
            //        *itLeft);
            else if (*itLeft < *itRight)
                ++itLeft;
            else
                ++itRight;
        }
        return product;
}


template <class T, class I>
void CSRUtility<T,I>::symmetricProduct(
    Ref<std::vector<I>> rTargetRowExtents,
    Ref<std::vector<I>> rTargetColumnIndices,
    Ref<std::vector<T>> rTargetEntries,
    CSRView<const T,const I> matrix,
    T scale) {
        // Step 1: compute the product's topology.
        CIE_BEGIN_EXCEPTION_TRACING
            std::vector<tsl::robin_set<I>> rows(matrix.rowCount());
            std::vector<std::mutex> mutexes(std::min<std::size_t>(
                matrix.rowCount(),
                _rMaybeThreads.has_value()
                    ? 1
                    : 1e2 * _rMaybeThreads.value().size()));

            const auto job = [&] (I iRow, Ref<std::vector<I>> rBuffer) -> void {
                const I iRowEntryBegin = matrix.rowExtents()[iRow];
                const I iRowEntryEnd   = matrix.rowExtents()[iRow + 1];
                for (I iEntry=iRowEntryBegin; iEntry<iRowEntryEnd; ++iEntry) {
                    const I row = matrix.columnIndices()[iEntry];
                    rBuffer.clear();
                    for (I jEntry=iRowEntryBegin; jEntry<iRowEntryEnd; ++jEntry) {
                        const I column = matrix.columnIndices()[jEntry];
                        const I iMutex = column % mutexes.size();
                        std::scoped_lock<std::mutex> lock(mutexes[iMutex]);
                        rows[column].insert(row);
                    } // for iEntry in range(iRowEntryBegin, iRowEntryEnd)
                    const I iMutex = row % mutexes.size();
                    std::scoped_lock<std::mutex> lock(mutexes[iMutex]);
                    rows[row].insert(
                        rBuffer.begin(),
                        rBuffer.end());
                } // for iEntry in range(iRowEntryBegin, iRowEntryEnd)
            }; // job

            // Execute the job that fills the hash map format of the matrix' graph.
            if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                mp::ParallelFor<I>(_rMaybeThreads.value()).firstPrivate(std::vector<I>()).execute(
                    matrix.rowCount(),
                    job);
            } else {
                std::vector<I> buffer;
                for (I iRow=0; iRow<matrix.rowCount(); ++iRow)
                    job(iRow, buffer);
            }

            // Build the matrix' CSR graph.
            rTargetRowExtents.resize(matrix.rowCount() + 1);
            rTargetRowExtents.front() = I(0);
            for (I iRow=0; iRow<static_cast<I>(rows.size()); ++iRow) {
                const auto& rRow = rows[iRow];
                rTargetRowExtents[iRow + 1] = rRow.size();
            } // for itRow
            std::inclusive_scan(
                rTargetRowExtents.begin(),
                rTargetRowExtents.end(),
                rTargetRowExtents.begin());
            rTargetColumnIndices.resize(rTargetRowExtents.back());
            rTargetEntries.resize(rTargetRowExtents.back());

            if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                    matrix.rowCount(),
                    [&] (const I iRow) -> void {
                        const auto& rRow = rows[iRow];
                        const Ptr<I> itTargetBegin = rTargetColumnIndices.data() + rTargetRowExtents[iRow];
                        const Ptr<I> itTargetEnd   = rTargetColumnIndices.data() + rTargetRowExtents[iRow + 1];
                        std::copy(
                            rRow.begin(),
                            rRow.end(),
                            itTargetBegin);
                        std::sort(
                            itTargetBegin,
                            itTargetEnd);});
            } /*if threads*/ else {
                for (I iRow=0; iRow<static_cast<I>(rows.size()); ++iRow) {
                    const auto& rRow = rows[iRow];
                    const Ptr<I> itTargetBegin = rTargetColumnIndices.data() + rTargetRowExtents[iRow];
                    const Ptr<I> itTargetEnd   = rTargetColumnIndices.data() + rTargetRowExtents[iRow + 1];
                    std::copy(
                        rRow.begin(),
                        rRow.end(),
                        itTargetBegin);
                    std::sort(
                        itTargetBegin,
                        itTargetEnd);
                }
            }
        CIE_END_EXCEPTION_TRACING

        CSRView<T,const I> target(
            matrix.columnCount(),
            rTargetRowExtents,
            rTargetColumnIndices,
            rTargetEntries);

        // Step 2: compute the product.
        CIE_BEGIN_EXCEPTION_TRACING
            const auto job = [&] (const I iEntry, auto&& op) -> void {
                // First of all, find which row and column the entry belongs to.
                const I iColumn = target.columnIndices()[iEntry];
                Ptr<const I> itRowExtent = std::upper_bound(
                    target.rowExtents().data(),
                    target.rowExtents().data() + target.rowExtents().size(),
                    iEntry) - 1;
                CIE_CHECK(
                    itRowExtent != target.rowExtents().data() + target.rowExtents().size(),
                    std::format(
                        "cannot find row that contains entry {}",
                        iEntry))
                const I iRow = std::distance<Ptr<const I>>(target.rowExtents().data(), itRowExtent);

                T product = 0;
                const I iRowEntryBegin    = matrix.rowExtents()[iRow];
                const I iRowEntryEnd      = matrix.rowExtents()[iRow + 1];
                if (iRow == iColumn) [[unlikely]] {
                    product = std::inner_product(
                        matrix.entries().data() + iRowEntryBegin,
                        matrix.entries().data() + iRowEntryEnd,
                        matrix.entries().data() + iRowEntryBegin,
                        T(0));
                } /*if iRow == iColumn*/ else [[likely]] {
                    product = this->sparseInnerProduct(iRow, iColumn, matrix);
                } /*if iRow != iColumn*/

                const T targetEntry = op(product);

                // Find the entry in the target matrix and set it.
                {
                    target.entries()[iEntry] = targetEntry;
                }
            }; // job

            if (_rMaybeThreads.has_value() && 1 < _rMaybeThreads.value().size()) {
                if (scale == static_cast<T>(1)) {
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        target.entries().size(),
                        [&job] (I iEntry) {job(iEntry, [] (T v) {return v;});});
                } else {
                    mp::ParallelFor<I>(_rMaybeThreads.value()).execute(
                        target.entries().size(),
                        [&job, scale] (I iEntry) {job(iEntry, [scale] (T v) {return scale * v;});});
                }
            } /*if threads*/ else {
                if (scale == static_cast<T>(1))
                    for (I iEntry=0; iEntry<static_cast<I>(target.entries().size()); ++iEntry)
                        job(iEntry, [] (T v) {return v;});
                else
                    for (I iEntry=0; iEntry<static_cast<I>(target.entries().size()); ++iEntry)
                        job(iEntry, [scale] (T v) {return scale * v;});
            } /*if not threads*/
        CIE_END_EXCEPTION_TRACING
}


#define CIE_DEFINE_CSR_UTILITY(T, I)        \
    template class CSRUtility<T,I>;

CIE_DEFINE_CSR_UTILITY(float, int)
CIE_DEFINE_CSR_UTILITY(float, std::size_t)
CIE_DEFINE_CSR_UTILITY(double, int)
CIE_DEFINE_CSR_UTILITY(double, std::size_t)

#undef CIE_DEFINE_CSR_UTILITY


} // namespace cie::linalg
