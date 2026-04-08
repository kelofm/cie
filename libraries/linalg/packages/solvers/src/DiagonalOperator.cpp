// --- Linalg Includes ---
#include "packages/solvers/inc/DiagonalOperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <optional>
#include <format>
#include <cassert>


namespace cie::linalg {


template <LinalgSpaceLike TS>
DiagonalOperator<TS>::DiagonalOperator()
requires std::is_default_constructible_v<typename TS::Vector>
    :   _pSpace(),
        _buffer(),
        _inverseDiagonal()
{}


template <LinalgSpaceLike TS>
DiagonalOperator<TS>::DiagonalOperator(
    typename TS::Vector&& rInverseDiagonal,
    std::shared_ptr<const TS> pSpace) noexcept
    :   _pSpace(pSpace),
        _buffer(),
        _inverseDiagonal(std::move(rInverseDiagonal)) {
    _buffer = pSpace->makeVector(pSpace->size(_inverseDiagonal));
}


template <LinalgSpaceLike TS>
void DiagonalOperator<TS>::product(
    typename TS::Value inScale,
    typename TS::ConstVectorView in,
    typename TS::Value outScale,
    typename TS::VectorView out) {
        using V = typename TS::Value;
        if (inScale == static_cast<V>(0)) {
            _pSpace->assign(out, in);
            _pSpace->scale(out, _inverseDiagonal, outScale);
        } else {
            _pSpace->assign(_buffer, in);
            _pSpace->scale(_buffer, _inverseDiagonal, outScale);
            if (inScale != static_cast<V>(1)) _pSpace->scale(out, inScale);
            _pSpace->add(out, _buffer, 1);
        }
}


template <
    class TV,
    class TI,
    class TMV>
DiagonalOperator<DefaultSpace<TV,tags::SMP>> makeDiagonalOperator(
    CSRView<const TMV,const TI> matrix,
    std::shared_ptr<const DefaultSpace<TV,tags::SMP>> pSpace) {
        assert(!matrix.rowExtents().empty());
        const TI rowCount = matrix.rowCount();
        auto inverseDiagonal = pSpace->makeVector(rowCount);

        for (TI iRow=0; iRow<rowCount; ++iRow) {
            const TI iEntryBegin = matrix.rowExtents()[iRow];
            const TI iEntryEnd = matrix.rowExtents()[iRow + 1];
            std::optional<TMV> maybeDiagonal;
            for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry)
                if (matrix.columnIndices()[iEntry] == iRow) [[unlikely]] {
                    maybeDiagonal = matrix.entries()[iEntry];
                    break;
                }
            if (maybeDiagonal.has_value()) [[likely]] {
                const TMV diagonal = *maybeDiagonal;
                CIE_CHECK(diagonal, std::format("row {} has an explicit zero on its diagonal", iRow))
                inverseDiagonal[iRow] = static_cast<TMV>(1) / diagonal;
            } else [[unlikely]] {
                CIE_THROW(Exception, std::format("row {} has no diagonal entry", iRow))
            }
        } // for iRow in range(iRow)

        return DiagonalOperator<DefaultSpace<TV,tags::SMP>>(
            std::move(inverseDiagonal),
            pSpace);
}


template class DiagonalOperator<DefaultSpace<float,tags::Serial>>;
template class DiagonalOperator<DefaultSpace<float,tags::SMP>>;
template class DiagonalOperator<DefaultSpace<double,tags::Serial>>;
template class DiagonalOperator<DefaultSpace<double,tags::SMP>>;


#define CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(TV, TI, TMV)                                   \
    template DiagonalOperator<DefaultSpace<TV,tags::SMP>> makeDiagonalOperator<TV,TI,TMV>(  \
        CSRView<const TMV,const TI>,                                                        \
        std::shared_ptr<const DefaultSpace<TV,tags::SMP>>);


CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(float, int, float)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(float, std::size_t, float)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(double, int, double)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(double, std::size_t, double)


#undef CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY


#ifdef CIE_ENABLE_SYCL


template <
    class TV,
    class TI,
    class TMV>
DiagonalOperator<SYCLSpace<TV>> makeDiagonalOperator(
    CSRView<const TMV,const TI> matrix,
    std::shared_ptr<const SYCLSpace<TV>> pSpace) {
        assert(!matrix.rowExtents().empty());

        // Make sure the matrix is on the device.
        Ptr<const TI> pRowExtentBegin = matrix.rowExtents().data();
        CIE_CHECK(
            sycl::get_pointer_type(pRowExtentBegin, pSpace->getQueue()->get_context()) == sycl::usm::alloc::device,
            "input matrix is not allocated on a SYCL device")

        const TI rowCount = matrix.rowCount();
        auto inverseDiagonal = pSpace->makeVector(rowCount);
        Ptr<TV> pInverseDiagonalBegin = inverseDiagonal.get();

        pSpace->getQueue()->parallel_for(
            sycl::range<1>(rowCount),
            [=] (sycl::item<1> it) -> void {
                const std::size_t iRow = it.get_linear_id();
                const std::size_t iEntryBegin = matrix.rowExtents()[iRow];
                const std::size_t iEntryEnd = matrix.rowExtents()[iRow + 1];

                if (iEntryEnd == iEntryBegin) {
                    pInverseDiagonalBegin[iRow] = std::numeric_limits<TV>::max();
                } else {
                    const Ptr<const TI> pColumnIndexBegin = matrix.columnIndices().data() + iEntryBegin;
                    const Ptr<const TI> pColumnIndexEnd = matrix.columnIndices().data() + iEntryEnd;
                    Ptr<const TI> pColumnIndex = std::upper_bound(
                        pColumnIndexBegin,
                        pColumnIndexEnd,
                        static_cast<TI>(iRow),
                        [] (TI iRow, TI iColumn) {return iRow < iColumn;});
                    if (pColumnIndex == pColumnIndexBegin) {
                        pInverseDiagonalBegin[iRow] = std::numeric_limits<TV>::max();
                    } else {
                        const std::size_t iEntry = std::distance(matrix.columnIndices().data(), pColumnIndex) - 1;
                        pInverseDiagonalBegin[iRow] = static_cast<TMV>(1) / matrix.entries()[iEntry];
                    }
                }

            }).wait_and_throw(); // parallel_for

        return DiagonalOperator<SYCLSpace<TV>>(
            std::move(inverseDiagonal),
            pSpace);
}


template class DiagonalOperator<SYCLSpace<float>>;
template class DiagonalOperator<SYCLSpace<double>>;


#define CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(TV, TI, TMV)                       \
    template DiagonalOperator<SYCLSpace<TV>> makeDiagonalOperator<TV,TI,TMV>(   \
        CSRView<const TMV,const TI>,                                            \
        std::shared_ptr<const SYCLSpace<TV>>);


CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(float, int, float)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(float, std::size_t, float)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(double, int, double)
CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY(double, std::size_t, double)


#undef CIE_DEFINE_DIAGONAL_OPERATOR_FACTORY


#endif


} // namespace cie::linalg
