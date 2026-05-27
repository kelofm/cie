// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- Linalg Includes ---
#include "packages/utilities/inc/CSRSum.hpp"

// --- STL Includes ---
#include <format>
#include <algorithm>


namespace cie::linalg {


template <class TValue, class TIndex>
void CSRSum(
    CSRView<TValue,TIndex> left,
    CSRView<const TValue,const TIndex> right,
    TValue scale,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads) {
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
            const auto makeJob = [left, right] (const auto& op) {
                return [left, right, op] (const TIndex iRow) mutable -> void {
                    const TIndex iLeftEntryBegin        = left.rowExtents()[iRow];
                    const TIndex iLeftEntryEnd          = left.rowExtents()[iRow + 1];
                    const Ptr<TIndex> itLeftColumnBegin = left.columnIndices().data() + iLeftEntryBegin;
                    const Ptr<TIndex> itLeftColumnEnd   = left.columnIndices().data() + iLeftEntryEnd;

                    const TIndex iRightEntryBegin = right.rowExtents()[iRow];
                    const TIndex iRightEntryEnd   = right.rowExtents()[iRow + 1];

                    Ptr<TIndex> itLeftColumn = itLeftColumnBegin;
                    for (TIndex iRightEntry=iRightEntryBegin; iRightEntry<iRightEntryEnd; ++iRightEntry) {
                        const TIndex iColumn = right.columnIndices()[iRightEntry];
                        itLeftColumn = std::lower_bound(
                            itLeftColumn,
                            itLeftColumnEnd,
                            iColumn);
                        CIE_CHECK(
                            itLeftColumn != itLeftColumnEnd && (*itLeftColumn) == iColumn,
                            std::format(
                                "left hand operand in matrix sum has no entry in row {} at column {}",
                                iRow, iColumn))
                        const TIndex iLeftEntry = static_cast<TIndex>(std::distance<Ptr<const TIndex>>(
                            left.columnIndices().data(),
                            itLeftColumn));
                        left.entries()[iLeftEntry] += op(right.entries()[iRightEntry]);
                    } // for iRightEntry in range(iRightEntryBegin, iRightEntryEnd)
                }; // job
            }; // makeJob

            if (rMaybeThreads.has_value() && 1 < rMaybeThreads.value().size()) {
                if (scale == static_cast<TValue>(1))
                    mp::ParallelFor<TIndex>(rMaybeThreads.value()).execute(
                        left.rowCount(),
                        makeJob([] (TValue r) -> TValue {return r;}));
                else
                    mp::ParallelFor<TIndex>(rMaybeThreads.value()).execute(
                        left.rowCount(),
                        makeJob([scale] (TValue r) -> TValue {return scale * r;}));
            } else {
                if (scale == static_cast<TValue>(1)) {
                    auto job = makeJob([] (TValue r) -> TValue {return r;});
                    for (TIndex iRow=0; iRow<left.rowCount(); ++iRow)
                        job(iRow);
                } else {
                    auto job = makeJob([scale] (TValue r) -> TValue {return scale * r;});
                    for (TIndex iRow=0; iRow<left.rowCount(); ++iRow)
                        job(iRow);
                }

            }
        CIE_END_EXCEPTION_TRACING
}


#define CIE_DEFINE_CSR_SUM(TV, TI)          \
    template void CSRSum<TV,TI>(            \
        CSRView<TV,TI>,                     \
        CSRView<const TV,const TI>,         \
        TV,                                 \
        OptionalRef<mp::ThreadPoolBase>);

CIE_DEFINE_CSR_SUM(float, int)
CIE_DEFINE_CSR_SUM(float, std::size_t)
CIE_DEFINE_CSR_SUM(double, int)
CIE_DEFINE_CSR_SUM(double, std::size_t)

#undef CIE_DEFINE_CSR_SUM


} // namespace cie::linalg
