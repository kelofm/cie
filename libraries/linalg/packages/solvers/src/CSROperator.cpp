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
    CSRView<const TMV,const TI> lhs,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads)
        :   _lhs(lhs),
            _maybeThreads(rMaybeThreads) {
    CIE_CHECK(
        0 < _lhs.rowExtents().size(),
        std::format(
            "expecting the row extents of a CSR matrix to have at least 1 entry, but it has none"))
    CIE_CHECK(
        _lhs.columnIndices().size() == _lhs.entries().size(),
        std::format(
            "number of column indices ({}) is inconsistent with the number of entries ({}) in the provided CSR matrix",
            _lhs.columnIndices().size(), _lhs.entries().size()))
}


template <class TI, class TV, class TMV>
void CSROperator<TI,TV,TMV>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        // Sanity checks.
        CIE_CHECK(
            _lhs.rowExtents().size() - 1 == out.size() && in.size() == static_cast<std::size_t>(_lhs.columnCount()),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] * [{}] = [{}]",
                _lhs.rowExtents().size() - 1, _lhs.columnCount(), in.size(), out.size()))

        const auto kernel = [this, in, out] (TI iRowBegin, TI iRowEnd, const auto& op) -> void {
            for (TI iRow=iRowBegin; iRow<iRowEnd; ++iRow) {
                const TI iEntryBegin = _lhs.rowExtents()[iRow];
                const TI iEntryEnd   = _lhs.rowExtents()[iRow + 1];
                TV contribution = static_cast<TV>(0);
                for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                    const TI iColumn = _lhs.columnIndices()[iEntry];
                    const TV entry   = _lhs.entries()[iEntry];
                    contribution += entry * in[iColumn];
                } // for iEntry in range(iEntryBegin, iEntryEnd)
                op(out[iRow], contribution);
            }
        }; // kernel

        const TI rowCount = _lhs.rowExtents().size() - 1;
        const auto job = [&kernel, rowCount, this] (const auto& op) -> void {
            if (_maybeThreads.has_value()) {
                mp::ParallelFor<TI>(_maybeThreads.value())
                    .execute(
                        mp::DynamicIndexPartitionFactory(
                            {0, static_cast<std::size_t>(rowCount), 1},
                            0x10 * _maybeThreads.value().size()),
                        [&kernel, &op] (
                            std::size_t iRowBegin,
                            std::size_t iRowEnd) -> void {
                                kernel(iRowBegin, iRowEnd, op);
                        });
                } else {
                    kernel(0, rowCount, op);
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


template class CSROperator<int,float,float>;
template class CSROperator<std::size_t,float,float>;
template class CSROperator<int,double,double>;
template class CSROperator<std::size_t,double,double>;


} // namespace cie::linalg
