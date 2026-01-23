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
    LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>,
    SpatialTransform::Derivative>;


static_assert(maths::StaticExpression<StiffnessIntegrand>);
static_assert(!StiffnessIntegrand::isBuffered);


class CellExtents {
public:
    using Value = std::tuple<
        VertexID,
        std::size_t,
        StiffnessIntegrand>;

    CellExtents()
        : _extents(SYCLSingleton::makeSharedAllocator<Value>())
    {
        this->clear();
    }

    std::size_t size() const noexcept {
        return std::get<std::size_t>(_extents.back());
    }

    bool empty() const noexcept {
        //return std::get<std::size_t>(_extents.back()) == 0ul;
        return _extents.size() == 1ul && std::get<std::size_t>(_extents.back()) == 0ul;
    }

    void push(VertexID cellID, RightRef<StiffnessIntegrand> rIntegrand) {
        _extents.push_back(std::make_tuple(
            cellID,
            std::get<std::size_t>(_extents.back()),
            std::move(rIntegrand)));
        std::swap(
            _extents[_extents.size() - 2],
            _extents[_extents.size() - 1]);
    }

    void extend(std::size_t count) noexcept {
        std::get<std::size_t>(_extents.back()) += count;
    }

    void clear() noexcept {
        _extents.resize(1);
        std::get<std::size_t>(_extents.back()) = 0ul;
        std::get<VertexID>(_extents.back()) = std::numeric_limits<unsigned>::max();
    }

    std::span<const Value> get() const noexcept {
        return _extents;
    }

private:
    DynamicSharedArray<Value> _extents;
}; // class CellExtents


struct QuadraturePointData {
    QuadraturePoint<Dimension,Scalar> point;
    std::array<Scalar,StiffnessIntegrand::size()> output;
};


} // namespace cie::fem


namespace cie::fem {


auto makeJob(std::span<QuadraturePointData> quadraturePointData,
             std::span<const CellExtents::Value> cellExtents,
             std::size_t integrandsPerItem = 1) {
    return [
        pQuadraturePointBegin   = quadraturePointData.data(),
        pCellExtentBegin        = cellExtents.data(),
        pCellExtentEnd          = cellExtents.data() + cellExtents.size(),
        quadraturePointCount    = quadraturePointData.size(),
        integrandsPerItem]
            (auto index) -> void {
                // Compute linear index.
                using TIndex = std::remove_cvref_t<decltype(index)>;
                std::size_t iItem = 0ul;
                if constexpr (concepts::UnsignedInteger<TIndex>) {
                    iItem = index;
                } else {
                    iItem = index.get_linear_id();
                }

                // Map to quadrature point range begin.
                const std::size_t iQuadraturePointBegin = iItem * integrandsPerItem;
                const std::size_t iQuadraturePointEnd = iQuadraturePointBegin + integrandsPerItem;

                for (std::size_t iQuadraturePoint=iQuadraturePointBegin; iQuadraturePoint<iQuadraturePointEnd; ++iQuadraturePoint) {
                    if (quadraturePointCount <= iQuadraturePoint) break;

                    Ref<const QuadraturePoint<Dimension,Scalar>> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint].point;

                    // Find which integrand the given quadrature point belongs to.
                    auto pExtent = std::upper_bound(
                        pCellExtentBegin,
                        pCellExtentEnd,
                        iQuadraturePoint,
                        [](std::size_t iQuadraturePoint, Ref<const CellExtents::Value> rExtent) -> bool {
                            return iQuadraturePoint < std::get<std::size_t>(rExtent);
                        });
                    if (iQuadraturePoint < std::get<std::size_t>(*pExtent)) --pExtent;

                    // Evaluate the integrand.
                    const StiffnessIntegrand integrand = std::get<StiffnessIntegrand>(*pExtent);
                    std::span<Scalar> output(
                        pQuadraturePointBegin[iQuadraturePoint].output.data(),
                        StiffnessIntegrand::size());
                    rQuadraturePoint.evaluate(integrand, output);
                }
            };
}


auto makeSingleJob(std::span<QuadraturePointData> quadraturePointData,
                   std::span<const CellExtents::Value> cellExtents) {
    return [
        pQuadraturePointBegin = quadraturePointData.data(),
        pCellExtentBegin = cellExtents.data(),
        pCellExtentEnd = cellExtents.data() + cellExtents.size()]
            (auto index) -> void {
                using TIndex = std::remove_cvref_t<decltype(index)>;
                std::size_t iQuadraturePoint = 0ul;
                if constexpr (concepts::UnsignedInteger<TIndex>) {
                    iQuadraturePoint = index;
                } else {
                    iQuadraturePoint = index.get_linear_id();
                }
                Ref<const QuadraturePoint<Dimension,Scalar>> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint].point;

                // Find which integrand the given quadrature point belongs to.
                auto pExtent = std::upper_bound(
                    pCellExtentBegin,
                    pCellExtentEnd,
                    iQuadraturePoint,
                    [](std::size_t iQuadraturePoint, Ref<const CellExtents::Value> rExtent) -> bool {
                        return iQuadraturePoint < std::get<std::size_t>(rExtent);
                    });
                if (iQuadraturePoint < std::get<std::size_t>(*pExtent)) --pExtent;

                // Evaluate the integrand.
                const StiffnessIntegrand integrand = std::get<StiffnessIntegrand>(*pExtent);
                std::span<Scalar> output(
                    pQuadraturePointBegin[iQuadraturePoint].output.data(),
                    StiffnessIntegrand::size());
                rQuadraturePoint.evaluate(integrand, output);
            };
}


class IntegrandRange {
public:
    IntegrandRange()
        : _quadraturePointData(SYCLSingleton::makeSharedAllocator<QuadraturePointData>()),
          _extents()
    {}

    void setQuadraturePointCount(std::size_t count) {
        _quadraturePointData.resize(count);
    }

    Ref<QuadraturePoint<Dimension,Scalar>> getQuadraturePoint(std::size_t iPoint) {
        return _quadraturePointData[iPoint].point;
    }

    void evaluate(Ref<const Assembler> rAssembler,
                  Ref<mp::ThreadPoolBase> rThreads,
                  bool useSYCL,
                  CSRWrapper lhs) {
        const auto cellExtents = _extents.get();

        if (cellExtents.empty()) return;
        auto logBlock = utils::LoggerSingleton::get().newBlock(
            "evaluate integrand at "
            + std::to_string(std::get<std::size_t>(cellExtents.back()))
            + " quadrature points");
        std::cout << "device memory footprint "
                  << std::get<std::size_t>(cellExtents.back()) * sizeof(QuadraturePointData)
                     + cellExtents.size() * sizeof(decltype(cellExtents)::value_type)
                  << "\n";

        #ifdef CIE_ENABLE_SYCL
            if (useSYCL) {
//                constexpr std::size_t integrandsPerItem = 0x10;
                SYCLSingleton::getQueue().parallel_for(
//                    sycl::range(_quadraturePointData.size() / integrandsPerItem + bool(_quadraturePointData.size() % integrandsPerItem)),
//                    makeJob(
//                        _quadraturePointData,
//                        cellExtents,
//                        integrandsPerItem)
                    sycl::range(_quadraturePointData.size()),
                    makeSingleJob(_quadraturePointData, cellExtents)
                ).wait_and_throw();
            } else {
                mp::ParallelFor<>(rThreads).operator()(
                    _quadraturePointData.size(),
                    makeSingleJob(_quadraturePointData, cellExtents));
            }
        #else
            mp::ParallelFor<>(rThreads).operator()(
                _quadraturePointData.size(),
                makeSingleJob(_quadraturePointData, cellExtents));
        #endif

        // Reduce the results and map them to the LHS matrix.
        constexpr unsigned integrandOutputSize = StiffnessIntegrand::size();
        for (std::size_t iCell=0ul; iCell<cellExtents.size()-1; ++iCell) {
            const VertexID cellID = std::get<VertexID>(cellExtents[iCell]);
            const std::size_t iBegin = std::get<std::size_t>(cellExtents[iCell]),
                              iEnd   = std::get<std::size_t>(cellExtents[iCell + 1]);
            const auto itTargetBegin = _quadraturePointData[iBegin].output.data();

            for (std::size_t iSample=iBegin+1; iSample<iEnd; ++iSample) {
                const auto itSourceBegin = _quadraturePointData[iSample].output.data();
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

    Ref<CellExtents> extents() noexcept {
        return _extents;
    }

    Ref<const CellExtents> extents() const noexcept {
        return _extents;
    }

private:
    DynamicSharedArray<QuadraturePointData> _quadraturePointData;

    CellExtents _extents;
}; // struct IntegrandDataRange


void integrateStiffness(Ref<const Mesh> rMesh,
                        Ref<const Assembler> rAssembler,
                        CSRWrapper lhs,
                        Ref<const utils::ArgParse::Results> rArguments,
                        Ref<mp::ThreadPoolBase> rThreads) {
    auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");
    const std::size_t quadratureBatchSize = rArguments.get<std::size_t>("integrand-batch-size");
    const bool useSYCL = rArguments.get<bool>("gpu");

    // Allocate buffers.
    IntegrandRange integrandRange;
    integrandRange.setQuadraturePointCount(quadratureBatchSize);

    // Perform integration and assembly.
    // 0) Generate quadrature points.
    // 1) Evaluate integrand at quadrature points.

    CachedQuadraturePointFactory<Dimension,Scalar> quadraturePointFactory;
    DynamicArray<QuadraturePoint<Dimension,Scalar>> quadraturePoints(quadratureBatchSize);

    for (const auto& rCell : rMesh.vertices()) {
        integrandRange.extents().push(
            rCell.data().id(),
            makeTransformedIntegrand(
                LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(
                    rCell.data().diffusivity(),
                    Ansatz::Derivative(rMesh.data().ansatzDerivative())),
                rCell.data().makeJacobian()));

        quadraturePointFactory = rMesh.data().makeQuadratureRule();
        while (true) {
            if (integrandRange.extents().empty()) {
                integrandRange.extents().push(
                    rCell.data().id(),
                    makeTransformedIntegrand(
                        LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(
                            rCell.data().diffusivity(),
                            Ansatz::Derivative(rMesh.data().ansatzDerivative())),
                        rCell.data().makeJacobian()));
            }

            // Fetch the range of remaining quadrature point slots.
            const std::size_t iIntegrandBegin = integrandRange.extents().size();
            std::span<QuadraturePoint<Dimension,Scalar>> quadraturePointSpan(
                quadraturePoints.data() + iIntegrandBegin,
                quadraturePoints.data() + quadraturePoints.size());

            // Generate new quadrature points.
            const std::size_t quadraturePointCount = quadraturePointFactory.generate(rCell.data(), quadraturePointSpan);
            const std::size_t iIntegrandEnd = iIntegrandBegin + quadraturePointCount;
            integrandRange.extents().extend(quadraturePointCount);

            for (std::size_t iIntegrand=iIntegrandBegin; iIntegrand<iIntegrandEnd; ++iIntegrand) {
                integrandRange.getQuadraturePoint(iIntegrand) = quadraturePoints[iIntegrand];
            } // for iIntegrand in range(iIntegrandBegin, iIntegrandEnd)

            if (quadraturePointCount == 0ul) {
                if (integrandRange.extents().size() == quadratureBatchSize) {
                    integrandRange.evaluate(
                        rAssembler,
                        rThreads,
                        useSYCL,
                        lhs);
                    integrandRange.extents().clear();
                } else break;
            }
        } // while true
    } // for rCell in rMesh.vertices()

    if (!integrandRange.extents().empty()) {
        integrandRange.setQuadraturePointCount(integrandRange.extents().size());
        integrandRange.evaluate(
            rAssembler,
            rThreads,
            useSYCL,
            lhs);
    }
}


} // namespace cie::fem
