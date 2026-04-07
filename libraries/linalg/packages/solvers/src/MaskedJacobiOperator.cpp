// --- Linalg Includes ---
#include "packages/solvers/inc/MaskedJacobiOperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <format>
#include <optional>



namespace cie::linalg {


template <class TI, class TV, class TMV, class TMI>
MaskedJacobiOperator<TI,TV,TMV,TMI>::MaskedJacobiOperator(
    TI columnCount,
    std::span<const TI> rowExtents,
    std::span<const TI> columnIndices,
    std::span<const TMV> entries,
    std::size_t iterations,
    TV relaxation,
    std::span<const TMI> mask,
    TMI threshold,
    std::shared_ptr<Space> pSpace)
        :   _columnCount(columnCount),
            _rowExtents(rowExtents),
            _columnIndices(columnIndices),
            _entries(entries),
            _mask(mask),
            _threshold(threshold),
            _iterations(iterations),
            _relaxation(relaxation),
            _previous(),
            _pSpace(pSpace) {
    CIE_CHECK(
        0 < _rowExtents.size(),
        std::format(
            "expecting the row extents of a CSR matrix to have at least 1 entry, but it has none"))
    CIE_CHECK(
        columnIndices.size() == entries.size(),
        std::format(
            "number of column indices ({}) is inconsistent with the number of entries ({}) in the provided CSR matrix",
            columnIndices.size(), entries.size()))
    _previous.resize(rowExtents.size() - 1);
}


template <class TI, class TV, class TMV, class TMI>
void MaskedJacobiOperator<TI,TV,TMV,TMI>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        // Sanity checks.
        CIE_CHECK(
            _rowExtents.size() - 1 == out.size() && in.size() == static_cast<std::size_t>(_columnCount),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] * [{}] = [{}]",
                _rowExtents.size() - 1, _columnCount, in.size(), out.size()))

        if (inScale != static_cast<TV>(0) && _iterations != 1)
            CIE_THROW(NotImplementedException, "")

        const auto kernel = [this, in, out] (
            TI iRowBegin,
            TI iRowEnd,
            const auto& op) -> void {
                for (TI iRow=iRowBegin; iRow<iRowEnd; ++iRow) {
                    if (_mask[iRow] < _threshold) {
                        TV contribution = static_cast<TV>(0);
                        const TI iEntryBegin = _rowExtents[iRow];
                        const TI iEntryEnd   = _rowExtents[iRow + 1];
                        std::optional<TMV> maybeDiagonal;
                        for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                            const TI iColumn = _columnIndices[iEntry];
                            const TV entry   = _entries[iEntry];
                            if (iRow != iColumn) [[likely]] {
                                if (_mask[iColumn] < _threshold)
                                    contribution += entry * _previous[iColumn];
                            } else {
                                maybeDiagonal = entry;
                            }
                        } // for iEntry in range(iEntryBegin, iEntryEnd)
                        op(out[iRow], _relaxation * (in[iRow] - contribution) / maybeDiagonal.value() + (static_cast<TV>(1) - _relaxation) * _previous[iRow]);
                    } else {
                        op(out[iRow], static_cast<TV>(0));
                    }
                }
        }; // kernel

        const TI rowCount = _rowExtents.size() - 1;
        const auto job = [&kernel, rowCount, out, this] (const auto& op) -> void {
            auto maybeThreads = _pSpace->getThreads();
            if (maybeThreads.has_value()) {
                mp::DynamicIndexPartitionFactory partitions(
                    {0, static_cast<std::size_t>(rowCount), 1},
                    0x10 * maybeThreads.value().size());
                auto loop = mp::ParallelFor<TI>(maybeThreads.value());
                for (std::size_t iIteration=0ul; iIteration<_iterations; ++iIteration) {
                    _pSpace->assign(_previous, out);
                    loop.execute(
                        partitions,
                        [&kernel, &op] (
                            std::size_t iRowBegin,
                            std::size_t iRowEnd) -> void {
                                kernel(iRowBegin, iRowEnd, op);});
                }
                } else {
                    for (std::size_t iIteration=0ul; iIteration<_iterations; ++iIteration) {
                        kernel(0, rowCount, op);
                        _pSpace->assign(_previous, out);
                    }
                }
        }; // job

        if (inScale == static_cast<TV>(1)) {
            if (outScale == static_cast<TV>(1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft += right;});
            } /*if outScale == 1*/ else if (outScale == static_cast<TV>(-1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft -= right;});
            } /*if outScale == -1*/ else {
                job([outScale] (Ref<TV> rLeft, TV right) -> void {rLeft += outScale * right;});
            }
        } /*if inScale == 1*/ else if (inScale == static_cast<TV>(0)) {
            if (outScale == static_cast<TV>(1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft = right;});
            } /*if outScale == 1*/ else if (outScale == static_cast<TV>(-1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft = -right;});
            } /*if outScale == -1*/ else {
                job([outScale] (Ref<TV> rLeft, TV right) -> void {rLeft = outScale * right;});
            }
        } /*if inScale == 0*/ else if (inScale == static_cast<TV>(-1)) {
            if (outScale == static_cast<TV>(1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft = right - rLeft;});
            } /*if outScale == 1*/ else if (outScale == static_cast<TV>(-1)) {
                job([] (Ref<TV> rLeft, TV right) -> void {rLeft = right - rLeft;});
            } /*if outScale == -1*/ else {
                job([outScale] (Ref<TV> rLeft, TV right) -> void {rLeft = outScale * right - rLeft;});
            }
        } /*if inScale == -1*/ else {
            if (outScale == static_cast<TV>(1)) {
                job([inScale] (Ref<TV> rLeft, TV right) -> void {rLeft = inScale * rLeft + right;});
            } /*if outScale == 1*/ else if (outScale == static_cast<TV>(-1)) {
                job([inScale] (Ref<TV> rLeft, TV right) -> void {rLeft = inScale * rLeft - right;});
            } /*if outScale == -1*/ else {
                job([inScale, outScale] (Ref<TV> rLeft, TV right) -> void {rLeft = inScale * rLeft + outScale * right;});
            }
        } // else
}


template class MaskedJacobiOperator<int,float,float,int>;
template class MaskedJacobiOperator<int,float,float,float>;
template class MaskedJacobiOperator<std::size_t,float,float,std::size_t>;
template class MaskedJacobiOperator<int,double,double,int>;
template class MaskedJacobiOperator<int,double,double,double>;
template class MaskedJacobiOperator<std::size_t,double,double,std::size_t>;


} // namespace cie::linalg
