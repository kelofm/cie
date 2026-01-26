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

    CellExtents(){
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

    std::span<Value> get() noexcept {
        return _extents;
    }

private:
    DynamicArray<Value> _extents;
}; // class CellExtents


//struct QuadraturePointData {
//     point;
//    //std::array<Scalar,StiffnessIntegrand::size()> output;
//};


} // namespace cie::fem


namespace cie::fem {


class IntegrandRange {
public:
    void setQuadraturePointCount(std::size_t count) {
        _quadraturePoints.resize(count);
    }

    Ref<QuadraturePoint<Dimension,Scalar>> getQuadraturePoint(std::size_t iPoint) {
        return _quadraturePoints[iPoint];
    }

    void evaluate(Ref<const Assembler> rAssembler,
                  Ref<mp::ThreadPoolBase> rThreads,
                  bool useSYCL,
                  CSRWrapper lhs) {
        if (_extents.get().empty()) return;
        auto logBlock = utils::LoggerSingleton::get().newBlock(
            "evaluate integrand at "
            + std::to_string(_quadraturePoints.size())
            + " quadrature points");

        if (useSYCL) {
            this->evaluateSYCL(rAssembler, lhs);
        } else {
            this->evaluateHost(rAssembler, rThreads, lhs);
        }
    }

    Ref<CellExtents> extents() noexcept {
        return _extents;
    }

    Ref<const CellExtents> extents() const noexcept {
        return _extents;
    }

private:
    auto makeJob(std::span<QuadraturePoint<Dimension,Scalar>> quadraturePoints,
                 std::span<const CellExtents::Value> cellExtents,
                 std::span<Scalar> output) {
            return [
                pQuadraturePointBegin   = quadraturePoints.data(),
                pCellExtentBegin        = cellExtents.data(),
                pCellExtentEnd          = cellExtents.data() + cellExtents.size(),
                pOutputBegin            = output.data(),
                quadraturePointCount    = quadraturePoints.size()]
                    (std::size_t iQuadraturePoint) -> void {
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
                        Ref<const QuadraturePoint<Dimension,Scalar>> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint];
                        const auto integrand = std::get<StiffnessIntegrand>(*pExtent);
                        std::span<Scalar> output(
                            pOutputBegin + iQuadraturePoint * StiffnessIntegrand::size(),
                            StiffnessIntegrand::size());
                        rQuadraturePoint.evaluate(integrand, output);
                    };
        }

    void evaluateHost(Ref<const Assembler> rAssembler,
                      Ref<mp::ThreadPoolBase> rThreads,
                      CSRWrapper lhs) {
        CIE_BEGIN_EXCEPTION_TRACING
        const auto cellExtents = _extents.get();
        Ptr<CellExtents::Value> pCellExtentBegin = cellExtents.data();
        Ptr<QuadraturePoint<Dimension,Scalar>> pQuadraturePointBegin = _quadraturePoints.data();
        DynamicArray<Scalar> output(_quadraturePoints.size() * StiffnessIntegrand::size());

        // Compute integrands.
        mp::ParallelFor<>(rThreads).operator()(
            _quadraturePoints.size(),
            makeJob(
                {pQuadraturePointBegin, _quadraturePoints.size()},
                {pCellExtentBegin, cellExtents.size()},
                output));

        // Reduce and assemble the results.
        constexpr unsigned integrandOutputSize = StiffnessIntegrand::size();
        for (std::size_t iCell=0ul; iCell<cellExtents.size()-1; ++iCell) {
            const VertexID cellID = std::get<VertexID>(cellExtents[iCell]);
            const std::size_t iBegin = std::get<std::size_t>(cellExtents[iCell]),
                                iEnd   = std::get<std::size_t>(cellExtents[iCell + 1]);
            const auto pTargetBegin = output.data() + iBegin * integrandOutputSize;

            for (std::size_t iSample=iBegin+1; iSample<iEnd; ++iSample) {
                const auto pSourceBegin = output.data() + iSample * integrandOutputSize;
                const auto pSourceEnd = pSourceBegin + integrandOutputSize;
                std::transform(
                    pSourceBegin,
                    pSourceEnd,
                    pTargetBegin,
                    pTargetBegin,
                    std::plus<Scalar>());
            } // for iSample in range(iBegin + 1, iEnd)

            // Map the reduced results to the LHS matrix.
            addLHSContribution(
                {pTargetBegin, integrandOutputSize},
                rAssembler[cellID],
                lhs);
        } // for iCell in range(extents.size())
        CIE_END_EXCEPTION_TRACING
    }

    #ifdef CIE_ENABLE_SYCL
        auto makeSYCLJob(std::span<QuadraturePoint<Dimension,Scalar>> quadraturePoints,
                         std::span<const CellExtents::Value> cellExtents,
                         std::span<Scalar> output,
                         std::size_t integrandsPerItem = 1) {
            return [
                pQuadraturePointBegin   = quadraturePoints.data(),
                pCellExtentBegin        = cellExtents.data(),
                pCellExtentEnd          = cellExtents.data() + cellExtents.size(),
                pOutputBegin            = output.data(),
                quadraturePointCount    = quadraturePoints.size(),
                integrandsPerItem]
                    (auto index) -> void {
                        // Compute linear index.
                        const std::size_t iItem = index.get_linear_id();

                        // Map to quadrature point range begin.
                        const std::size_t iQuadraturePointBegin = iItem * integrandsPerItem;
                        const std::size_t iQuadraturePointEnd = iQuadraturePointBegin + integrandsPerItem;
                        std::size_t iLastQuadraturePoint = std::numeric_limits<std::size_t>::max();
                        StiffnessIntegrand integrand;
                        std::array<Scalar,StiffnessIntegrand::size()> result;

                        for (std::size_t iQuadraturePoint=iQuadraturePointBegin; iQuadraturePoint<iQuadraturePointEnd; ++iQuadraturePoint) {
                            if (quadraturePointCount <= iQuadraturePoint) break;

                            Ref<const QuadraturePoint<Dimension,Scalar>> rQuadraturePoint = pQuadraturePointBegin[iQuadraturePoint];

                            // Find which integrand the given quadrature point belongs to.
                            auto pExtent = std::upper_bound(
                                pCellExtentBegin,
                                pCellExtentEnd,
                                iQuadraturePoint,
                                [](std::size_t iQuadraturePoint, Ref<const CellExtents::Value> rExtent) -> bool {
                                    return iQuadraturePoint < std::get<std::size_t>(rExtent);
                                });
                            if (iQuadraturePoint < std::get<std::size_t>(*pExtent)) --pExtent;

                            if (iLastQuadraturePoint != iQuadraturePoint)
                                integrand = std::get<StiffnessIntegrand>(*pExtent);

                            // Evaluate the integrand.
                            rQuadraturePoint.evaluate(integrand, result);

                            std::span<Scalar> output(
                                pOutputBegin + std::distance(pCellExtentBegin, pExtent) * StiffnessIntegrand::size(),
                                StiffnessIntegrand::size());
                            for (std::size_t iComponent=0ul; iComponent<StiffnessIntegrand::size(); ++iComponent) {
                                sycl::atomic_ref<
                                    Scalar,
                                    sycl::memory_order::relaxed,
                                    sycl::memory_scope::device
                                >(output[iComponent]) += result[iComponent];
                            }

                            iLastQuadraturePoint = iQuadraturePoint;
                        }
                    };
        }
    #endif

    void evaluateSYCL(Ref<const Assembler> rAssembler,
                      CSRWrapper lhs) {
        #ifdef CIE_ENABLE_SYCL
            auto cellExtents = _extents.get();
            std::cout << "device memory footprint "
                  << cellExtents.size() * sizeof(CellExtents::Value)
                     + _quadraturePoints.size() * sizeof(QuadraturePoint<Dimension,Scalar>)
                     + _quadraturePoints.size() * StiffnessIntegrand::size() * sizeof(Scalar)
                  << "\n";

            // Allocate device memory.
            Ref<sycl::queue> rQueue = SYCLSingleton::getQueue();
            Ptr<CellExtents::Value> pCellExtentBegin = sycl::malloc_device<CellExtents::Value>(
                    cellExtents.size(),
                    rQueue);
            Ptr<QuadraturePoint<Dimension,Scalar>> pQuadraturePointBegin = sycl::malloc_device<QuadraturePoint<Dimension,Scalar>>(
                    _quadraturePoints.size(),
                    rQueue);
            Ptr<Scalar> pOutput = sycl::malloc_device<Scalar>(
                    (cellExtents.size() - 1) * StiffnessIntegrand::size(),
                    rQueue);

            // Copy input data to the device.
            const std::vector<sycl::event> copyEvents {
                rQueue.submit([cellExtents, pCellExtentBegin](sycl::handler& rHandler) {
                    rHandler.copy(
                        cellExtents.data(),
                        pCellExtentBegin,
                        cellExtents.size());}),
                rQueue.submit([this, pQuadraturePointBegin](sycl::handler& rHandler) {
                    rHandler.copy(
                        _quadraturePoints.data(),
                        pQuadraturePointBegin,
                        _quadraturePoints.size());})};

            // Perform calculation on the device.
            constexpr std::size_t integrandsPerItem = 0x8;
            sycl::event computeEvent = rQueue.parallel_for(
                sycl::range(_quadraturePoints.size() / integrandsPerItem + bool(_quadraturePoints.size() % integrandsPerItem)),
                copyEvents,
                makeSYCLJob(
                    {pQuadraturePointBegin, _quadraturePoints.size()},
                    {pCellExtentBegin, cellExtents.size()},
                    {pOutput, (cellExtents.size() - 1) * StiffnessIntegrand::size()},
                    integrandsPerItem)
            );
            computeEvent.wait_and_throw();

            // Fetch results from the device.
            DynamicArray<Scalar> output((cellExtents.size() - 1) * StiffnessIntegrand::size());
            auto returnEvent = rQueue.submit([pOutput, &output, this](sycl::handler& rHandler){
                rHandler.copy(
                    pOutput,
                    output.data(),
                    (_extents.get().size() - 1) * StiffnessIntegrand::size());});

            // Reduce and assemble the results.
            constexpr unsigned integrandOutputSize = StiffnessIntegrand::size();
            for (std::size_t iCell=0ul; iCell<cellExtents.size()-1; ++iCell) {
                const VertexID cellID = std::get<VertexID>(cellExtents[iCell]);
                const auto pTargetBegin = output.data() + iCell * integrandOutputSize;

                // Map the reduced results to the LHS matrix.
                addLHSContribution(
                    {pTargetBegin, integrandOutputSize},
                    rAssembler[cellID],
                    lhs);
            } // for iCell in range(extents.size())

            // Release device memory.
            sycl::free(pCellExtentBegin, rQueue);
            sycl::free(pQuadraturePointBegin, rQueue);
            sycl::free(pOutput, rQueue);
        #else
            CIE_THROW(Exception, "Requested SYCL evaluation, but CiE was compiled without SYCL support.")
        #endif
    }

    DynamicArray<QuadraturePoint<Dimension,Scalar>> _quadraturePoints;

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
