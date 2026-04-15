// --- Internal Includes ---
#include "embeddedPoisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/ScaledMultiMaterialIntegrand.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"
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
    std::span<const CellData> contiguousCellData,
    Ref<const Assembler> rAssembler,
    linalg::CSRView<Scalar,int> lhs,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads) {
        auto logBlock = utils::LoggerSingleton::get().newBlock("integrate stiffness matrix");
        const std::size_t quadratureBatchSize = rArguments.get<std::size_t>("integrand-batch-size");

        using IntegrandBase = LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>;
        using EmbeddedIntegrand = ScaledMultiMaterialIntegrand<IntegrandBase,Scalar,2>;
        using Integrand = TransformedIntegrand<
            EmbeddedIntegrand,
            SpatialTransform::Inverse::Derivative>;
        static_assert(maths::StaticExpression<Integrand>);

        std::vector<std::unique_ptr<IntegrandProcessor<
            Dimension,
            Integrand,
            Scalar
        >>> integrandProcessors;
        std::vector<std::size_t> partitions;
        partitions.push_back(0ul);

        CIE_BEGIN_EXCEPTION_TRACING
        #ifdef CIE_ENABLE_SYCL
            if (rArguments.get<bool>("sycl")) {
                auto logBlock = utils::LoggerSingleton::get().newBlock("discover SYCL devices");
                std::vector<sycl::device> devices;
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
                        Integrand,
                        Scalar>>(std::make_shared<sycl::queue>(device)));
                } // for device in devices
                partitions.back() = contiguousCellData.size();
            } else {
                partitions.push_back(contiguousCellData.size());
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
            }
        #else
            partitions.push_back(contiguousCellData.size());
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
        const auto quadratureRuleFactory = [&rMesh] (Ref<const Mesh::Vertex::Data>) {
            return rMesh.data().makeQuadratureRule();};

        const auto integrandFactory = [&rMesh] (Ref<const Mesh::Vertex::Data> rCell) {
            const std::array<std::pair<Scalar,Scalar>,2> materialMap {
                std::pair<Scalar,Scalar> {0, std::numeric_limits<Scalar>::epsilon()},
                std::pair<Scalar,Scalar> {1, 1}};
            return Integrand(
                EmbeddedIntegrand(
                    IntegrandBase(
                        rCell.diffusivity(),
                        Ansatz::Derivative(rMesh.data().ansatzDerivative(rCell.ansatzID()))),
                    materialMap),
                rCell.makeJacobianInverse());};

        const auto integralSink = [&lhs, &rAssembler, &rThreads] (std::span<const VertexID> cellIDs, std::span<const Scalar> results) {
            mp::ParallelFor<std::size_t>(rThreads).execute(
                cellIDs.size(),
                [&lhs, &rAssembler, cellIDs, results] (std::size_t iCell) {
                    rAssembler.addContribution<tags::SMP,int>(
                        std::span<const Scalar>(results.data() + iCell * Integrand::size(), Integrand::size()),
                        cellIDs[iCell],
                        lhs.rowExtents(),
                        lhs.columnIndices(),
                        lhs.entries());
                });};

        IntegrandProcessor<Dimension,Integrand>::Properties executionProperties {
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
