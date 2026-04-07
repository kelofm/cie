// --- Linalg Includes ---
#include "packages/solvers/inc/ApproximateSubstitution.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <format>


namespace cie::linalg {


template <class TV, class TI, class TMV>
ApproximateSubstitution<TV,TI,TMV>::ApproximateSubstitution(
    CSRView<const TMV,const TI> lhs,
    std::size_t iterations,
    std::shared_ptr<Space> pSpace)
        :   _lhs(lhs),
            _pSpace(pSpace),
            _iterations(iterations),
            _previous(lhs.rowCount())
{}


template <class TV, class TI, class TMV>
void ApproximateSubstitution<TV,TI,TMV>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        // Sanity checks.
        CIE_CHECK(
            _lhs.rowCount() == out.size() && in.size() == static_cast<std::size_t>(_lhs.columnCount()),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] @ [{}] = [{}]",
                _lhs.rowCount(), _lhs.columnCount(), in.size(), out.size()))

        if (inScale != static_cast<TV>(0) && _iterations != 1)
            CIE_THROW(NotImplementedException, "")

        const auto kernel = [this, in, out] (TI iRow, const auto& op) -> void {
            const TI iEntryBegin = _lhs.rowExtents()[iRow];
            const TI iEntryEnd   = _lhs.rowExtents()[iRow + 1];
            std::optional<TMV> maybeDiagonal;
            TV contribution = static_cast<TV>(0);
            for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                const TI iColumn = _lhs.columnIndices()[iEntry];
                const TV entry   = _lhs.entries()[iEntry];
                if (iRow != iColumn) [[likely]]
                    contribution += entry * _previous[iColumn];
                else
                    maybeDiagonal = entry;
            } // for iEntry in range(iEntryBegin, iEntryEnd)
            const TV diagonal = static_cast<TV>(maybeDiagonal.value());
            op(out[iRow], _previous[iRow] + (in[iRow] - contribution) / diagonal);
        }; // kernel

        const auto job = [&kernel, out, this] (const auto& op) -> void {
            auto maybeThreads = _pSpace->getThreads();
            if (maybeThreads.has_value()) {
                auto loop = mp::ParallelFor<TI>(maybeThreads.value());
                for (std::size_t iIteration=0ul; iIteration<_iterations; ++iIteration) {
                    _pSpace->assign(_previous, out);
                    loop.execute(
                        _lhs.rowCount(),
                        [&kernel, &op] (std::size_t iRow) -> void {kernel(iRow, op);});
                }
                } else {
                    for (std::size_t iIteration=0ul; iIteration<_iterations; ++iIteration) {
                        for (TI iRow=0; iRow<_lhs.rowCount(); ++iRow) kernel(iRow, op);
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


} // namespace cie::linalg
