#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"
#include "packages/utilities/inc/IntegrandProcessor.hpp"

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


void integrateStiffness(Ref<const Mesh> rMesh,
                        Ref<const Assembler> rAssembler,
                        CSRWrapper lhs,
                        Ref<const utils::ArgParse::Results> rArguments,
                        Ref<mp::ThreadPoolBase> rThreads) {
    auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");
    const std::size_t quadratureBatchSize = rArguments.get<std::size_t>("integrand-batch-size");

    using StiffnessIntegrand = TransformedIntegrand<
        LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>,
        SpatialTransform::Inverse::Derivative>;
    static_assert(maths::StaticExpression<StiffnessIntegrand>);
    static_assert(!StiffnessIntegrand::isBuffered);

    auto pIntegrandProcessor = std::make_unique<IntegrandProcessor<
        Dimension,
        StiffnessIntegrand>>();

    CIE_BEGIN_EXCEPTION_TRACING
    #ifdef CIE_ENABLE_SYCL
        if (rArguments.get<bool>("gpu")) {
            pIntegrandProcessor = std::make_unique<SYCLIntegrandProcessor<
                Dimension,
                StiffnessIntegrand>>(std::make_shared<sycl::queue>());
        } else if (rThreads.size() != 1) {
            pIntegrandProcessor = std::make_unique<ParallelIntegrandProcessor<
                Dimension,
                StiffnessIntegrand>>(rThreads);
        }
    #else
        if (rThreads.size() != 1) {
            pIntegrandProcessor = std::make_unique<ParallelIntegrandProcessor<
                Dimension,
                StiffnessIntegrand>>(rThreads);
        }
    #endif
    CIE_END_EXCEPTION_TRACING

    CIE_BEGIN_EXCEPTION_TRACING
    const auto quadratureRuleFactory = [] (Ref<const Mesh> rMesh, Ref<const Mesh::Vertex::Data>) {
        return rMesh.data().makeQuadratureRule();};

    const auto integrandFactory = [] (Ref<const Mesh> rMesh, Ref<const CellData> rCell) {
        return StiffnessIntegrand(
            LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(
                rCell.diffusivity(),
                Ansatz::Derivative(rMesh.data().ansatzDerivative())),
            rCell.makeJacobianInverse());};

    const auto integralSink = [&lhs, &rAssembler, &rThreads] (std::span<const VertexID> cellIDs, std::span<const Scalar> results) {
        mp::ParallelFor<std::size_t>(rThreads).operator()(
            cellIDs.size(),
            [&lhs, &rAssembler, cellIDs, results] (std::size_t iCell) {
                rAssembler.addContribution(
                    std::span<const Scalar>(results.data() + iCell * StiffnessIntegrand::size(), StiffnessIntegrand::size()),
                    cellIDs[iCell],
                    lhs.rowExtents,
                    lhs.columnIndices,
                    lhs.entries);
            });};

    IntegrandProcessor<Dimension,StiffnessIntegrand>::Properties executionProperties {
        .integrandBatchSize = quadratureBatchSize,
        .integrandsPerItem = {}};

    pIntegrandProcessor->process(
        rMesh,
        quadratureRuleFactory,
        integrandFactory,
        integralSink,
        executionProperties);
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
