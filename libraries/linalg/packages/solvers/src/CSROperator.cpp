// --- Linalg Includes ---
#include "packages/solvers/inc/CSROperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <format>


namespace cie::linalg {


template <class TI, class TV, class TMV>
CSROperator<TI,TV,TMV>::CSROperator(
    TI columnCount,
    std::span<const TI> rowExtents,
    std::span<const TI> columnIndices,
    std::span<const TMV> entries,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads)
        :
        _columnCount(columnCount),
        _rowExtents(rowExtents),
        _columnIndices(columnIndices),
        _entries(entries),
        _maybeThreads(rMaybeThreads) {
    CIE_CHECK(
        0 < _rowExtents.size(),
        std::format(
            "expecting the row extents of a CSR matrix to have at least 1 entry, but it has none"))
    CIE_CHECK(
        columnIndices.size() == entries.size(),
        std::format(
            "number of column indices ({}) is inconsistent with the number of entries ({}) in the provided CSR matrix",
            columnIndices.size(), entries.size()))
}


template <class TI, class TV, class TMV>
void CSROperator<TI,TV,TMV>::product(
    typename Space::ConstVectorView in,
    typename Space::Value scale,
    typename Space::VectorView out) const {
        // Sanity checks.
        CIE_CHECK(
            _rowExtents.size() - 1 == out.size() && in.size() == static_cast<std::size_t>(_columnCount),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] * [{}] = [{}]",
                _rowExtents.size() - 1, _columnCount, in.size(), out.size()))

        const auto kernel = [this, in, out] (TI iRow, const auto& op) -> void {
            const TI iEntryBegin = _rowExtents[iRow];
            const TI iEntryEnd   = _rowExtents[iRow + 1];
            TV contribution = static_cast<TV>(0);
            for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                const TI iColumn = _columnIndices[iEntry];
                const TV entry   = _entries[iEntry];
                contribution += entry * in[iColumn];
            } // for iEntry in range(iEntryBegin, iEntryEnd)
            op(out[iRow], contribution);
        }; // kernel

        const TI rowCount = _rowExtents.size() - 1;

        if (scale == static_cast<TV>(1)) {
            if (_maybeThreads.has_value()) {
                mp::ParallelFor<TI>(_maybeThreads.value())
                    .firstPrivate([] (Ref<TV> rLeft, TV right) -> void {std::atomic_ref<TV>(rLeft) += right;})
                    .execute(
                        _rowExtents.size() - 1,
                        kernel);
            } else {
                const auto op = [] (Ref<TV> rLeft, TV right) -> void {rLeft += right;};
                for (TI iRow=0; iRow<rowCount; ++iRow)
                    kernel(iRow, op);
            }
        } /*if scale == 1*/ else if (scale == static_cast<TV>(-1)) {
            if (_maybeThreads.has_value()) {
                mp::ParallelFor<TI>(_maybeThreads.value())
                    .firstPrivate([] (Ref<TV> rLeft, TV right) -> void {std::atomic_ref<TV>(rLeft) -= right;})
                    .execute(
                        _rowExtents.size() - 1,
                        kernel);
            } else {
                const auto op = [] (Ref<TV> rLeft, TV right) -> void {rLeft -= right;};
                for (TI iRow=0; iRow<rowCount; ++iRow)
                    kernel(iRow, op);
            }
        } /*if scale == -1*/ else {
            if (_maybeThreads.has_value()) {
                const auto op = [scale] (Ref<TV> rLeft, TV right) -> void {std::atomic_ref<TV>(rLeft) += scale * right;};
                mp::ParallelFor<TI>(_maybeThreads.value())
                    .execute(
                        _rowExtents.size() - 1,
                        [&op, &kernel] (std::size_t iRow) -> void {kernel(iRow, op);});
            } else {
                const auto op = [scale] (Ref<TV> rLeft, TV right) -> void {rLeft += scale * right;};
                for (TI iRow=0; iRow<rowCount; ++iRow)
                    kernel(iRow, op);
            }
        }
}


template class CSROperator<int,float,float>;
template class CSROperator<std::size_t,float,float>;

template class CSROperator<int,double,double>;
template class CSROperator<std::size_t,double,double>;


} // namespace cie::linalg
