// --- Linalg Includes ---
#include "packages/solvers/inc/JacobiOperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <optional>
#include <format>
#include <cassert>


namespace cie::linalg {


template <LinalgSpaceLike TS>
JacobiOperator<TS>::JacobiOperator()
requires std::is_default_constructible_v<typename TS::Vector>
    :   _pSpace(),
        _inverseDiagonal()
{}


template <LinalgSpaceLike TS>
JacobiOperator<TS>::JacobiOperator(
    typename TS::Vector&& rInverseDiagonal,
    std::shared_ptr<const TS> pSpace) noexcept
    :   _pSpace(pSpace),
        _inverseDiagonal(std::move(rInverseDiagonal))
{}


template <LinalgSpaceLike TS>
void JacobiOperator<TS>::product(
    typename TS::ConstVectorView in,
    typename TS::Value scale,
    typename TS::VectorView out) const {
        _pSpace->assign(out, in);
        _pSpace->scale(
            out,
            _pSpace->view(_inverseDiagonal),
            scale);
}


template <
    class TV,
    class TI,
    class TMV>
JacobiOperator<DefaultSpace<TV,tags::SMP>> makeJacobiOperator(
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

        return JacobiOperator<DefaultSpace<TV,tags::SMP>>(
            std::move(inverseDiagonal),
            pSpace);
}


template class JacobiOperator<DefaultSpace<float,tags::Serial>>;
template class JacobiOperator<DefaultSpace<float,tags::SMP>>;
template class JacobiOperator<DefaultSpace<double,tags::Serial>>;
template class JacobiOperator<DefaultSpace<double,tags::SMP>>;

#define CIE_DEFINE_JACOBI_OPERATOR_FACTORY(TV, TI, TMV)                                                     \
    template JacobiOperator<DefaultSpace<TV,tags::SMP>> makeJacobiOperator<TV,TI,TMV>(                      \
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
