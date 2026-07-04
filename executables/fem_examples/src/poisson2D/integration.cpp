// --- Internal Includes ---
#include "poisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/ScaledMultiMaterialIntegrand.hpp"
#include "packages/utilities/inc/IntegrandProcessor.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <ranges>


namespace cie::fem {


void integrateStiffness(
    Ref<const Mesh> rMesh,
    std::span<const Cell> cells,
    Ref<const Assembler> rAssembler,
    linalg::CSRView<Scalar,int> lhs,
    Ref<const cie::io::JSONObject> rConfiguration,
    Ref<mp::ThreadPoolBase> rThreads) {
        auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");
        const std::size_t quadratureBatchSize = rConfiguration["batch-size"].as<std::size_t>();

        using IntegrandBase = LinearIsotropicStiffnessIntegrand<Ansatz::Derivative,SpatialTransform>;
        using Integrand = ScaledMultiMaterialIntegrand<IntegrandBase,Scalar,2>;
        static_assert(maths::StaticExpression<Integrand>);

        std::vector<std::unique_ptr<IntegrandProcessor<
            Dimension,
            Integrand,
            Scalar
        >>> integrandProcessors;
        std::vector<std::size_t> partitions;
        partitions.push_back(0ul);

        const std::string integrandProcessorDeviceName = rConfiguration["device"].as<std::string>();

        CIE_BEGIN_EXCEPTION_TRACING
            #ifdef CIE_ENABLE_SYCL
                if (integrandProcessorDeviceName == "sycl") {
                    auto logBlock = utils::LoggerSingleton::get().newBlock("discover SYCL devices");
                    std::vector<sycl::device> devices;
                    for (auto device : sycl::device::get_devices(sycl::info::device_type::gpu)) {
                        std::cout << device.get_info<sycl::info::device::name>() << std::endl;
                        devices.push_back(device);
                    } // for device
                    for (auto device : devices) {
                        partitions.push_back(std::min(
                            partitions.back() + cells.size() / devices.size(),
                            cells.size()));
                        integrandProcessors.emplace_back(std::make_unique<SYCLIntegrandProcessor<
                            Dimension,
                            Integrand,
                            Scalar>>(std::make_shared<sycl::queue>(device)));
                    } // for device in devices
                    partitions.back() = cells.size();
                } else if (integrandProcessorDeviceName == "host") {
                    partitions.push_back(cells.size());
                    if (rThreads.size() == 1) {
                        integrandProcessors.emplace_back(std::make_unique<IntegrandProcessor<
                                Dimension,
                                Integrand,
                                Scalar>>());
                    } else {
                        integrandProcessors.emplace_back(std::make_unique<ParallelIntegrandProcessor<
                            Dimension,
                            Integrand,
                            Scalar>>(rThreads));
                    }
                } else CIE_THROW(Exception, std::format(
                    "unsupported device \"{}\" for integration",
                    integrandProcessorDeviceName))
            #else
                CIE_CHECK(
                    integrandProcessorDeviceName == "host",
                    std::format(
                        "unsupported device \"{}\" for integration",
                        integrandProcessorDeviceName))
                partitions.push_back(cells.size());
                if (rThreads.size() == 1) {
                    integrandProcessors.emplace_back(std::make_unique<IntegrandProcessor<
                            Dimension,
                            Integrand,
                            Scalar
                        >>());
                } else {
                    integrandProcessors.emplace_back(std::make_unique<ParallelIntegrandProcessor<
                        Dimension,
                        Integrand,
                        Scalar>>(rThreads));
                }
            #endif
        CIE_END_EXCEPTION_TRACING

        CIE_BEGIN_EXCEPTION_TRACING
            const auto quadratureRuleFactory = [&rMesh] (Ref<const Cell> rCell) {
                return rMesh.data().makeQuadratureRule(rCell);};

            const auto integrandFactory = [&rMesh] (Ref<const Cell> rCell) -> Integrand {
                CIE_CHECK(rMesh.data().domainMap().size() == 2, "")
                std::span<const std::pair<MeshData::DomainData,Scalar>,2> domainData(
                    rMesh.data().domainMap().data(),
                    2);
                return Integrand(
                        IntegrandBase(
                            rCell.diffusivity(),
                            rMesh.data().ansatzDerivative(rCell.ansatzID()),
                            rCell.makeJacobian(),
                            rCell.makeJacobianInverse()),
                        domainData);};

            const auto integralSink = [&lhs, &rAssembler, &rThreads] (std::span<const VertexID> cellIDs, std::span<const Scalar> results) {
                mp::ParallelFor<std::size_t>(rThreads).execute(
                    cellIDs.size(),
                    [&lhs, &rAssembler, cellIDs, results] (std::size_t iCell) {
                        rAssembler.addContribution<tags::SMP,int>(
                            std::span<const Scalar>(
                                results.data() + iCell * Integrand::size(),
                                Integrand::size()),
                            cellIDs[iCell],
                            lhs.rowExtents(),
                            lhs.columnIndices(),
                            lhs.entries());
                    });};

            IntegrandProcessor<Dimension,Integrand,Scalar>::Properties executionProperties {
                .integrandBatchSize = quadratureBatchSize,
                .verbosity = rConfiguration["verbosity"].as<int>()};

            {
                std::vector<std::thread> jobs;
                jobs.reserve(integrandProcessors.size());
                for (std::size_t iPartition=0ul; iPartition<integrandProcessors.size(); ++iPartition) {
                    jobs.emplace_back([&, iPartition] () {
                        integrandProcessors[iPartition]->process(
                            std::span<const Cell> (
                                cells.begin() + partitions[iPartition],
                                cells.begin() + partitions[iPartition + 1]),
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
