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
    _buffer = pSpace->makeVector(pSpace->size(rInverseDiagonal));
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
    std::span<const TI> rowExtents,
    std::span<const TI> columnIndices,
    std::span<const TMV> entries,
    std::shared_ptr<const DefaultSpace<TV,tags::SMP>> pSpace) {
        assert(!rowExtents.empty());
        const TI rowCount = rowExtents.size() - 1;
        auto inverseDiagonal = pSpace->makeVector(rowCount);

        for (TI iRow=0; iRow<rowCount; ++iRow) {
            const TI iEntryBegin = rowExtents[iRow];
            const TI iEntryEnd = rowExtents[iRow + 1];
            std::optional<TMV> maybeDiagonal;
            for (TI iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry)
                if (columnIndices[iEntry] == iRow) [[unlikely]] {
                    maybeDiagonal = entries[iEntry];
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


#define CIE_DEFINE_JACOBI_OPERATOR_FACTORY(TV, TI, TMV)                                                     \
    template DiagonalOperator<DefaultSpace<TV,tags::SMP>> makeDiagonalOperator<TV,TI,TMV>(                      \
        std::span<const TI>,                                                                                \
        std::span<const TI>,                                                                                \
        std::span<const TMV>,                                                                               \
        std::shared_ptr<const DefaultSpace<TV,tags::SMP>>);


CIE_DEFINE_JACOBI_OPERATOR_FACTORY(float, int, float)
CIE_DEFINE_JACOBI_OPERATOR_FACTORY(float, std::size_t, float)
CIE_DEFINE_JACOBI_OPERATOR_FACTORY(double, int, double)
CIE_DEFINE_JACOBI_OPERATOR_FACTORY(double, std::size_t, double)


#undef CIE_DEFINE_JACOBI_OPERATOR_FACTORY

} // namespace cie::linalg
