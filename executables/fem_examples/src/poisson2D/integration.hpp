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
                        std::span<const CellData> contiguousCellData,
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

    std::vector<std::unique_ptr<IntegrandProcessor<
        Dimension,
        StiffnessIntegrand
    >>> integrandProcessors;
    std::vector<std::size_t> partitions {0ul};

    CIE_BEGIN_EXCEPTION_TRACING
    #ifdef CIE_ENABLE_SYCL
        if (rArguments.get<bool>("sycl")) {
            auto logBlock = utils::LoggerSingleton::get().newBlock("discover SYCL devices");
            std::vector<sycl::device> devices;
            for (auto device : sycl::device::get_devices(sycl::info::device_type::cpu)) {
                std::cout << device.get_info<sycl::info::device::name>() << std::endl;
                devices.push_back(device);
                break;
            } // for device
            for (auto device : sycl::device::get_devices(sycl::info::device_type::gpu)) {
                std::cout << device.get_info<sycl::info::device::name>() << std::endl;
                devices.push_back(device);
            } // for device
            for (auto device : devices) {
                partitions.push_back(std::min(
                    partitions.back() + contiguousCellData.size() / devices.size(),
                    contiguousCellData.size()));
                integrandProcessors.emplace_back(std::make_unique<SYCLIntegrandProcessor<
                    Dimension,
                    StiffnessIntegrand>>(std::make_shared<sycl::queue>(device)));
            } // for device in devices
            partitions.back() = contiguousCellData.size();
        } else {
            partitions.push_back(contiguousCellData.size());
            if (rThreads.size() == 1) {
                integrandProcessors.emplace_back(std::make_unique<IntegrandProcessor<
                        Dimension,
                        StiffnessIntegrand>>());
            } else {
                integrandProcessors.emplace_back(std::make_unique<ParallelIntegrandProcessor<
                    Dimension,
                    StiffnessIntegrand>>(rThreads));
            }
        }
    #else
        if (rThreads.size() == 1) {
            integrandProcessors.emplace_back(std::make_unique<IntegrandProcessor<
                    Dimension,
                    StiffnessIntegrand
                >>());
        } else {
            integrandProcessors.emplace_back(std::make_unique<ParallelIntegrandProcessor<
                Dimension,
                StiffnessIntegrand>>(rThreads));
        }
    #endif
    CIE_END_EXCEPTION_TRACING

    CIE_BEGIN_EXCEPTION_TRACING
    const auto quadratureRuleFactory = [&rMesh] (Ref<const Mesh::Vertex::Data>) {
        return rMesh.data().makeQuadratureRule();};

    const auto integrandFactory = [&rMesh] (Ref<const Mesh::Vertex::Data> rCell) {
        return StiffnessIntegrand(
            LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(
                rCell.diffusivity(),
                Ansatz::Derivative(rMesh.data().ansatzDerivative())),
            rCell.makeJacobianInverse());};

    std::mutex integralSinkMutex;
    const auto integralSink = [&lhs, &rAssembler, &rThreads, &integralSinkMutex] (std::span<const VertexID> cellIDs, std::span<const Scalar> results) {
        std::scoped_lock<std::mutex> lock(integralSinkMutex);
        mp::ParallelFor<std::size_t>(rThreads).operator()(
            cellIDs.size(),
            [&lhs, &rAssembler, cellIDs, results] (std::size_t iCell) {
                rAssembler.addContribution<tags::SMP>(
                    std::span<const Scalar>(results.data() + iCell * StiffnessIntegrand::size(), StiffnessIntegrand::size()),
                    cellIDs[iCell],
                    lhs.rowExtents,
                    lhs.columnIndices,
                    lhs.entries);
            });};

    IntegrandProcessor<Dimension,StiffnessIntegrand>::Properties executionProperties {
        .integrandBatchSize = quadratureBatchSize,
        .integrandsPerItem = {},
        .verbosity = 3};

    {
        std::vector<std::thread> jobs;
        jobs.reserve(integrandProcessors.size());
        for (std::size_t iPartition=0ul; iPartition<integrandProcessors.size(); ++iPartition) {
            jobs.emplace_back([&, iPartition] () {
                integrandProcessors[iPartition]->process(
                    contiguousCellData.begin() + partitions[iPartition],
                    contiguousCellData.begin() + partitions[iPartition + 1],
                    quadratureRuleFactory,
                    integrandFactory,
                    integralSink,
                    executionProperties);
                });
        }
        for (Ref<std::thread> rJob : jobs) rJob.join();
        jobs.clear();
    }
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
