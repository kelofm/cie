#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <ranges>



namespace cie::fem {


struct CSRWrapper {
    const int rowCount, columnCount;
    std::span<const int> rowExtents, columnIndices;
    std::span<Scalar> entries;
}; // struct CSRWrapper


template <class TDofMap>
void addLHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        CSRWrapper lhs) {
    const unsigned localSystemSize = dofMap.size();
    for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
        for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
            const auto iRowBegin = lhs.rowExtents[dofMap[iLocalRow]];
            const auto iRowEnd = lhs.rowExtents[dofMap[iLocalRow] + 1];
            const auto itColumnIndex = std::lower_bound(lhs.columnIndices.begin() + iRowBegin,
                                                        lhs.columnIndices.begin() + iRowEnd,
                                                        dofMap[iLocalColumn]);
            CIE_OUT_OF_RANGE_CHECK(itColumnIndex != lhs.columnIndices.begin() + iRowEnd
                                && *itColumnIndex == dofMap[iLocalColumn]);
            const auto iEntry = std::distance(lhs.columnIndices.begin(),
                                            itColumnIndex);
            lhs.entries[iEntry] += contribution[iLocalRow * localSystemSize + iLocalColumn];
        } // for iLocalColumn in range(ansatzBuffer.size)
    } // for iLocalRow in range(ansatzBuffer.size)
}


template <class TDofMap>
void addRHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        std::span<Scalar> rhs) {
    for (unsigned iComponent=0u; iComponent<contribution.size(); ++iComponent) {
        const auto iRow = dofMap[iComponent];
        rhs[iRow] += contribution[iComponent];
    }
}


using StiffnessIntegrand = TransformedIntegrand<
    LinearIsotropicStiffnessIntegrand<AnsatzDerivative>,
    SpatialTransform::Derivative
>;


class CellExtents {
public:
    CellExtents()
        : _extents()
    {
        this->clear();
    }

    std::size_t size() const noexcept {
        return _extents.back().second;
    }

    bool empty() const noexcept {
        return _extents.back().second == 0ul;
    }

    void push(VertexID cellID) {
        _extents.push_back(std::make_pair(
            cellID,
            _extents.back().second
        ));
        std::swap(
            _extents[_extents.size() - 2],
            _extents[_extents.size() - 1]);
    }

    void extend(std::size_t count) noexcept {
        _extents.back().second += count;
    }

    void clear() noexcept {
        _extents.resize(1);
        _extents.back() = std::make_pair(
            VertexID(std::numeric_limits<unsigned>::max()),
            0ul);
    }

    std::span<const std::pair<
        VertexID,
        std::size_t
    >> get() const noexcept {
        return _extents;
    }

private:
    DynamicArray<std::pair<
        VertexID,
        std::size_t
    >> _extents;
}; // class CellExtents


class IntegrandRange {
public:
    IntegrandRange()
        : _quadraturePoints(SYCLSingleton::makeSharedAllocator<QuadraturePoint<Dimension,Scalar>>()),
          _integrands(SYCLSingleton::makeSharedAllocator<StiffnessIntegrand>()),
          _integrandOutput(SYCLSingleton::makeSharedAllocator<Scalar>()),
          _integrandBuffer(SYCLSingleton::makeSharedAllocator<Scalar>())
    {}

    void setQuadraturePointCount(std::size_t count,
                                 std::size_t integrandOutputSize,
                                 std::size_t integrandBufferSize) {
        _quadraturePoints.resize(count);
        _integrands.resize(count);
        _integrandOutput.resize(count * integrandOutputSize);
        _integrandBuffer.resize(count * integrandBufferSize);
    }

    std::span<QuadraturePoint<Dimension,Scalar>> getQuadraturePoints(std::size_t iBegin) {
        CIE_OUT_OF_RANGE_CHECK(iBegin < _quadraturePoints.size())
        return {
            _quadraturePoints.data() + iBegin,
            _quadraturePoints.data() + _quadraturePoints.size()};
    }

    void makeIntegrand(Ref<const CellData> rCell,
                       Ref<const Mesh> rMesh,
                       std::size_t iIntegrand,
                       std::size_t integrandBufferSize) {
        _integrands[iIntegrand] = makeTransformedIntegrand(
            LinearIsotropicStiffnessIntegrand<AnsatzDerivative>(
                rCell.diffusivity(),
                AnsatzDerivative(rMesh.data().ansatzDerivatives()[rCell.ansatzID()]),
                {_integrandBuffer.data() + iIntegrand * integrandBufferSize, integrandBufferSize}),
            rCell.makeJacobian());
    }

    void evaluate(std::span<const std::pair<VertexID,std::size_t>> cellExtents,
                  Ref<const Assembler> rAssembler,
                  [[maybe_unused]] OptionalRef<mp::ThreadPoolBase> rMaybeThreads,
                  CSRWrapper lhs) {
        if (cellExtents.empty()) return;
        auto logBlock = utils::LoggerSingleton::get().newBlock(
            "evaluate integrand at "
            + std::to_string(cellExtents.back().second)
            + " quadrature points");

        const std::size_t integrandOutputSize = _integrands.empty()
            ? 0ul
            : _integrands.front().size();

        // Evaluate integrands at quadrature points.
        const auto job = [
            pQuadraturePointBegin = _quadraturePoints.data(),
            pIntegrandBegin = _integrands.data(),
            pIntegrandOutputBegin = _integrandOutput.data(),
            integrandOutputSize]
                (auto index) -> void {
                    std::size_t iQuadraturePoint = 0ul;
                    if constexpr (concepts::UnsignedInteger<decltype(index)>) {
                        iQuadraturePoint = index;
                    } else {
                        iQuadraturePoint = index.get_linear_id();
                    }
                    Ref<const QuadraturePoint<Dimension,Scalar>> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint];
                    Ref<const StiffnessIntegrand> rIntegrand = pIntegrandBegin[iQuadraturePoint];
                    std::span<Scalar> output(
                        pIntegrandOutputBegin + iQuadraturePoint * integrandOutputSize,
                        integrandOutputSize);
                    rQuadraturePoint.evaluate(rIntegrand, output);
                };

        if (rMaybeThreads.has_value()) {
            mp::ParallelFor<>(rMaybeThreads.value())(
                _quadraturePoints.size(),
                job);
        } else {
            #ifdef CIE_ENABLE_SYCL
                SYCLSingleton::getQueue().parallel_for(
                    sycl::range(_quadraturePoints.size()),
                    job
                ).wait();
            #else
                std::ranges::for_each(
                    std::views::iota(0ul) | std::views::take(_quadraturePoints.size()),
                    job);
            #endif
        }

        // Reduce the results and map them to the LHS matrix.
        for (std::size_t iCell=0ul; iCell<cellExtents.size()-1; ++iCell) {
            const VertexID cellID = cellExtents[iCell].first;
            const std::size_t iBegin = cellExtents[iCell].second,
                              iEnd   = cellExtents[iCell + 1].second;
            const auto itTargetBegin = _integrandOutput.data() + iBegin * integrandOutputSize;

            for (std::size_t iSample=iBegin+1; iSample<iEnd; ++iSample) {
                const auto itSourceBegin = _integrandOutput.data() + iSample * integrandOutputSize;
                const auto itSourceEnd = itSourceBegin + integrandOutputSize;
                std::transform(
                    itSourceBegin,
                    itSourceEnd,
                    itTargetBegin,
                    itTargetBegin,
                    std::plus<Scalar>());
            } // for iSample in range(iBegin + 1, iEnd)

            // Map the reduced results to the LHS matrix.
            addLHSContribution(
                {itTargetBegin, integrandOutputSize},
                rAssembler[cellID],
                lhs);
        } // for iCell in range(extents.size())
    }

private:
    DynamicSharedArray<QuadraturePoint<Dimension,Scalar>> _quadraturePoints;
    DynamicSharedArray<StiffnessIntegrand> _integrands;
    DynamicSharedArray<Scalar> _integrandOutput;
    DynamicSharedArray<Scalar> _integrandBuffer;
}; // struct IntegrandDataRange


void integrateStiffness(Ref<const Mesh> rMesh,
                        Ref<const Assembler> rAssembler,
                        CSRWrapper lhs,
                        std::size_t quadratureBatchSize,
                        OptionalRef<mp::ThreadPoolBase> rMaybeThreads) {
    auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");

    // Query integrand memory requirements.
    std::size_t integrandOutputSize, integrandBufferSize;
    if (!rMesh.vertices().empty()) {
        Ref<const CellData> rCell = rMesh.vertices().front().data();
        const StiffnessIntegrand integrand = makeTransformedIntegrand(
            LinearIsotropicStiffnessIntegrand<AnsatzDerivative>(
                rCell.diffusivity(),
                AnsatzDerivative(rMesh.data().ansatzDerivatives()[rCell.ansatzID()])),
            rCell.makeJacobian());
        integrandOutputSize = integrand.size();
        integrandBufferSize = integrand.getMinBufferSize();
    } else return;

    // Allocate buffers.
    IntegrandRange integrandRange;
    integrandRange.setQuadraturePointCount(
        quadratureBatchSize,
        integrandOutputSize,
        integrandBufferSize);

    // Perform integration and assembly.
    // 0) Generate quadrature points.
    // 1) Evaluate integrand at quadrature points.

    CachedQuadraturePointFactory<Dimension,Scalar> quadraturePointFactory;
    CellExtents extents;

    for (const auto& rCell : rMesh.vertices()) {
        extents.push(rCell.id());
        quadraturePointFactory = rMesh.data().makeQuadratureRule(rCell.data().ansatzID());
        while (true) {
            if (extents.empty())
                extents.push(rCell.id());

            // Fetch the range of remaining quadrature point slots.
            auto quadraturePoints = integrandRange.getQuadraturePoints(extents.size());
            const std::size_t iIntegrandBegin = extents.size();

            // Generate new quadrature points.
            const std::size_t quadraturePointCount = quadraturePointFactory.generate(rCell.data(), quadraturePoints);
            const std::size_t iIntegrandEnd = iIntegrandBegin + quadraturePointCount;
            extents.extend(quadraturePointCount);

            for (std::size_t iIntegrand=iIntegrandBegin; iIntegrand<iIntegrandEnd; ++iIntegrand) {
                integrandRange.makeIntegrand(
                    rCell.data(),
                    rMesh,
                    iIntegrand,
                    integrandBufferSize);
            } // for iIntegrand in range(iIntegrandBegin, iIntegrandEnd)

            if (quadraturePointCount == 0ul) {
                if (extents.size() == quadratureBatchSize) {
                    integrandRange.evaluate(
                        extents.get(),
                        rAssembler,
                        rMaybeThreads,
                        lhs);
                    extents.clear();
                } else break;
            }
        } // while true
    } // for rCell in rMesh.vertices()

    if (!extents.empty()) {
        integrandRange.setQuadraturePointCount(
            extents.size(),
            integrandOutputSize,
            integrandBufferSize);
        integrandRange.evaluate(
            extents.get(),
            rAssembler,
            rMaybeThreads,
            lhs);
    }
}


} // namespace cie::fem
