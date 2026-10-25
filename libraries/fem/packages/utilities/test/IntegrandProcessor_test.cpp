// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/ScaledMultiMaterialIntegrand.hpp"
#include "packages/utilities/inc/IntegrandProcessor.hpp"

#ifdef CIE_ENABLE_SYCL
    // --- Utility Includes ---
    #include "packages/concurrency/inc/sycl.hpp"
    #include "packages/stl_extension/inc/LinearAllocator.hpp"
    #include "packages/stl_extension/inc/TrackedAllocator.hpp"

    // --- STL Includes ---
    #include <memory>
    #include <sstream>
    #include <string>
#endif // CIE_ENABLE_SYCL


namespace cie::fem {


class ActiveTestCell
    : public CellBase<1,double,maths::IdentityTransform<double,1>> {
public:
    using Base = CellBase<1,double,maths::IdentityTransform<double,1>>;
    using Context = HyperSlab<1,std::uint16_t>;

    ActiveTestCell(VertexID id, Context activeContext)
        : Base(
            id,
            AnsatzID(0u),
            OrientedAxes<1>(),
            maths::IdentityTransform<double,1>()),
          _activeContext(activeContext)
    {}

    Context activeContext() const noexcept {
        return _activeContext;
    }

private:
    Context _activeContext;
}; // class ActiveTestCell


#ifdef CIE_ENABLE_SYCL
    class ActiveTestCell2D
        : public CellBase<2,double,maths::IdentityTransform<double,2>> {
    public:
        using Base = CellBase<2,double,maths::IdentityTransform<double,2>>;
        using Context = HyperSlab<2,std::uint16_t>;

        ActiveTestCell2D(VertexID id, Context activeContext)
            : Base(
                id,
                AnsatzID(0u),
                OrientedAxes<2>(),
                maths::IdentityTransform<double,2>()),
            _activeContext(activeContext)
        {}

        Context activeContext() const noexcept {
            return _activeContext;
        }

    private:
        Context _activeContext;
    }; // class ActiveTestCell2D
#endif // CIE_ENABLE_SYCL


CIE_TEST_CASE("IntegrandProcessor tiled slice assembly", "[utilities][integrands]") {
    CIE_TEST_CASE_INIT("IntegrandProcessor tiled slice assembly")

    using Basis = maths::Polynomial<double>;
    using Integrand = maths::AnsatzSpaceView<Basis,1>;
    using Point = QuadraturePoint<1,double>;
    using Processor = IntegrandProcessor<1,Integrand>;

    const std::array<Basis,3> basisFunctions {
        Basis({1.0}),
        Basis({0.5,0.5}),
        Basis({2.0})};
    const Integrand fullIntegrand(basisFunctions);

    std::array<double,1> position {0.0};
    const std::span<const double,1> positionSpan(position);
    const std::array<Point,2> quadraturePoints {
        Point(positionSpan,0.5),
        Point(positionSpan,0.5)};

    ActiveTestCell::Context activeContext;
    activeContext.data()[0] = 1;
    activeContext.data()[1] = 3;
    const std::array<ActiveTestCell,1> cells {
        ActiveTestCell(VertexID(7u),activeContext)};

    const auto quadratureRuleFactory = [&quadraturePoints] (Ref<const ActiveTestCell>) {
        return CachedQuadraturePointFactory<1,double>(quadraturePoints);
    };
    const auto integrandFactory = [&basisFunctions] (Ref<const ActiveTestCell>) {
        return Integrand(basisFunctions);
    };

    std::vector<std::vector<double>> batchResults;
    const auto integralSink = [&batchResults] (
        std::span<const VertexID> cellIDs,
        std::span<const double> values) {
            CIE_TEST_CHECK(cellIDs.size() == 1ul);
            CIE_TEST_CHECK(values.size() == 3ul);
            batchResults.emplace_back(values.begin(), values.end());
    };

    Processor processor;
    Processor::Properties properties {
        .integrandBatchSize = 1ul,
        .verbosity = 0u};
    properties.maxTileSizes[0] = 1ul;
    processor.process(
        std::span<const ActiveTestCell>(cells),
        quadratureRuleFactory,
        integrandFactory,
        integralSink,
        properties);

    // Two points with a batch size of one produce two superposable partial outputs.
    CIE_TEST_REQUIRE(batchResults.size() == 2ul);
    std::array<double,3> fullValues;
    std::vector<double> buffer(fullIntegrand.bufferSize());
    fullIntegrand.evaluate(positionSpan, fullValues, buffer);
    std::array<double,3> accumulated{};
    for (const auto& rBatch : batchResults) {
        CIE_TEST_CHECK(rBatch[0] == Approx(0.0)); // inactive function
        CIE_TEST_CHECK(rBatch[1] == Approx(0.5 * fullValues[1]));
        CIE_TEST_CHECK(rBatch[2] == Approx(0.5 * fullValues[2]));
        for (std::size_t iComponent=0ul; iComponent<accumulated.size(); ++iComponent)
            accumulated[iComponent] += rBatch[iComponent];
    }
    CIE_TEST_CHECK(accumulated[0] == Approx(0.0));
    CIE_TEST_CHECK(accumulated[1] == Approx(fullValues[1]));
    CIE_TEST_CHECK(accumulated[2] == Approx(fullValues[2]));

    // A two-context matrix integrand must visit row/column tile pairs, including
    // off-diagonal pairs, and scatter them into the dense cell matrix.
    using Ansatz = maths::AnsatzSpace<Basis,1>;
    using Transform = maths::IdentityTransform<double,1>;
    using LHSIntegrand = DirichletPenaltyLHSIntegrand<Ansatz,Transform,ActiveTestCell>;
    const Ansatz ansatzSpace(Ansatz::AnsatzSet(basisFunctions.begin(), basisFunctions.end()));
    const auto lhsFactory = [&] (Ref<const ActiveTestCell>) {
        const Transform transform;
        return LHSIntegrand(
            10.0,
            ansatzSpace,
            transform,
            transform.makeInverse(),
            transform.makeDerivative());
    };
    std::vector<std::vector<double>> matrixBatches;
    const auto matrixSink = [&matrixBatches] (
        std::span<const VertexID>,
        std::span<const double> values) {
            CIE_TEST_CHECK(values.size() == 9ul);
            matrixBatches.emplace_back(values.begin(), values.end());
    };
    const std::array<std::uint16_t,1> tileShape {1};
    auto pMatrixTiling = std::make_unique<typename IntegrandProcessor<1,LHSIntegrand>::ContextTiling>(tileShape);
    IntegrandProcessor<1,LHSIntegrand> matrixProcessor(std::move(pMatrixTiling));
    const IntegrandProcessor<1,LHSIntegrand>::Properties matrixProperties {
        .integrandBatchSize = 1ul,
        .verbosity = 0u};
    matrixProcessor.process(
        std::span<const ActiveTestCell>(cells),
        quadratureRuleFactory,
        lhsFactory,
        matrixSink,
        matrixProperties);

    CIE_TEST_REQUIRE(matrixBatches.size() == 2ul);
    const std::array<double,9> expectedMatrixContribution {
        0.0, 0.0, 0.0,
        0.0, 2.5, 10.0,
        0.0, 10.0, 40.0};
    std::array<double,9> accumulatedMatrix{};
    for (const auto& rBatch : matrixBatches) {
        for (std::size_t iComponent=0ul; iComponent<rBatch.size(); ++iComponent) {
            CIE_TEST_CHECK(rBatch[iComponent] == Approx(0.5 * expectedMatrixContribution[iComponent]));
            accumulatedMatrix[iComponent] += rBatch[iComponent];
        }
    }
    for (std::size_t iComponent=0ul; iComponent<accumulatedMatrix.size(); ++iComponent)
        CIE_TEST_CHECK(accumulatedMatrix[iComponent] == Approx(expectedMatrixContribution[iComponent]));

    using StaticBasis = maths::Polynomial<double,1>;
    using StaticIntegrand = maths::AnsatzSpaceView<StaticBasis,1,3>;
    const std::array<StaticBasis,3> staticBasisFunctions {
        StaticBasis(StaticBasis::Coefficients{1.0,0.0}),
        StaticBasis(StaticBasis::Coefficients{0.5,0.5}),
        StaticBasis(StaticBasis::Coefficients{2.0,0.0})};
    const auto staticIntegrandFactory = [&staticBasisFunctions] (Ref<const ActiveTestCell>) {
        return StaticIntegrand(std::span<const StaticBasis,3>(staticBasisFunctions));
    };
    std::vector<std::vector<double>> parallelBatches;
    const auto parallelSink = [&parallelBatches] (
        std::span<const VertexID>,
        std::span<const double> values) {
            parallelBatches.emplace_back(values.begin(), values.end());
    };
    const std::array<ActiveTestCell,2> parallelCells {
        cells.front(),
        ActiveTestCell(VertexID(8u),activeContext)};
    mp::ThreadPoolBase threads(2);
    ParallelIntegrandProcessor<1,StaticIntegrand> parallelProcessor(threads);
    const IntegrandProcessor<1,StaticIntegrand>::Properties parallelProperties {
        .integrandBatchSize = 1ul,
        .verbosity = 0u};
    parallelProcessor.process(
        std::span<const ActiveTestCell>(parallelCells),
        quadratureRuleFactory,
        staticIntegrandFactory,
        parallelSink,
        parallelProperties);
    CIE_TEST_REQUIRE(parallelBatches.size() == 4ul);
    for (std::size_t iBatch=0ul; iBatch<parallelBatches.size(); ++iBatch)
        CIE_TEST_REQUIRE(parallelBatches[iBatch] == batchResults[iBatch % batchResults.size()]);
}


#ifdef CIE_ENABLE_SYCL
    CIE_TEST_CASE("SYCLIntegrandProcessor dynamically sized integrands", "[utilities][integrands][sycl]") {
        CIE_TEST_CASE_INIT("SYCLIntegrandProcessor dynamically sized integrands")

        using Scalar = double;
        constexpr unsigned Dimension = 2u;

        // A dynamically sized integrand: the coefficients of its polynomial basis
        // functions and the basis sets themselves are heap allocated.
        using Basis = maths::Polynomial<Scalar>::Rebind<TrackedStandardAllocator,Scalar>;
        using Ansatz = maths::AnsatzSpace<Basis,Dimension>::Rebind<TrackedStandardAllocator,Scalar>;
        using Transform = maths::IdentityTransform<Scalar,Dimension>;
        using Base = LinearIsotropicStiffnessIntegrand<Ansatz::Derivative,Transform>;
        using Integrand = ScaledMultiMaterialIntegrand<Base,Scalar,1>;
        static_assert(!maths::StaticExpression<Integrand>);
        static_assert(!std::is_trivially_copyable_v<Integrand>);
        static_assert(cie::BinarySerializable<Integrand>);

        // Allocator statistics measuring the dynamic footprint of each integrand.
        AllocatorStats stats;
        TrackedStandardAllocator<Scalar> allocator(
            std::allocator<Scalar>(),
            stats);

        // Define a bilinear ansatz space over [-1,1]^2 (4 ansatz functions).
        const auto makeAnsatzDerivatives = [&] {
            typename Ansatz::AnsatzSet set({
                Basis(Basis::Coefficients({0.5, 0.5}, allocator)),
                Basis(Basis::Coefficients({0.5, -0.5}, allocator))},
                allocator);
            return Ansatz(std::move(set)).makeDerivative();
        };

        // Cells activating the complete ansatz space.
        ActiveTestCell2D::Context activeContext;
        activeContext.data()[0] = 0;
        activeContext.data()[1] = 2;
        activeContext.data()[2] = 0;
        activeContext.data()[3] = 2;
        const std::array<ActiveTestCell2D,2> cells {
            ActiveTestCell2D(VertexID(7u),activeContext),
            ActiveTestCell2D(VertexID(8u),activeContext)};

        // Four Gauss-Legendre points on [-1,1]^2.
        constexpr Scalar coordinate = 0.5773502691896257;
        const std::array<QuadraturePoint<Dimension,Scalar,Scalar>,4> quadraturePoints {
            QuadraturePoint<Dimension,Scalar,Scalar>(
                std::span<const Scalar,Dimension>(std::array<Scalar,Dimension> {-coordinate,-coordinate}.data(), Dimension),
                1.0,
                0.0),
            QuadraturePoint<Dimension,Scalar,Scalar>(
                std::span<const Scalar,Dimension>(std::array<Scalar,Dimension> { coordinate,-coordinate}.data(), Dimension),
                1.0,
                0.0),
            QuadraturePoint<Dimension,Scalar,Scalar>(
                std::span<const Scalar,Dimension>(std::array<Scalar,Dimension> {-coordinate, coordinate}.data(), Dimension),
                1.0,
                0.0),
            QuadraturePoint<Dimension,Scalar,Scalar>(
                std::span<const Scalar,Dimension>(std::array<Scalar,Dimension> { coordinate, coordinate}.data(), Dimension),
                1.0,
                0.0)};
        const auto quadratureRuleFactory = [&quadraturePoints] (Ref<const ActiveTestCell2D>) {
            return CachedQuadraturePointFactory<Dimension,Scalar,Scalar>(quadraturePoints);
        };

        // The integrand factory constructs integrands with the tracked allocator,
        // which SYCLIntegrandProcessor measures through the shared statistics.
        const std::array<std::pair<Scalar,Scalar>,1> materialMap {{{0.0,1.0}}};
        const auto integrandFactory = [&] (Ref<const ActiveTestCell2D>) -> Integrand {
            return Integrand(
                Base(
                    10.0,
                    makeAnsatzDerivatives(),
                    Transform().makeDerivative(),
                    Transform().makeInverse().makeDerivative()),
                std::span<const std::pair<Scalar,Scalar>,1>(materialMap));
        };

        // Verify the same serialized representation and allocator replacement on
        // the host first, using the device integrand type and a linear arena.
        {
            using DeviceIntegrand = typename Integrand::template Rebind<LinearAllocator,std::byte>;
            const Integrand source = integrandFactory(cells.front());
            std::ostringstream serializedStream;
            cie::io::BinarySerializer::serialize(serializedStream, source);
            const std::string serializedString = serializedStream.str();
            const std::span<const std::byte> serializedBytes(
                reinterpret_cast<const std::byte*>(serializedString.data()),
                serializedString.size());

            std::vector<std::byte> pool(4096ul);
            LinearArena arena { .buffer = pool, .state = 0ul };
            LinearAllocator<std::byte> linearAllocator(arena);
            DeviceIntegrand deserialized;
            cie::io::ByteDeserializer deserializer(serializedBytes);
            cie::io::BinarySerializer::deserialize(
                deserializer,
                deserialized,
                linearAllocator);

            std::vector<Scalar> sourceResult(source.size());
            std::vector<Scalar> deserializedResult(deserialized.size());
            std::vector<Scalar> sourceBuffer(source.bufferSize());
            std::vector<Scalar> deserializedBuffer(deserialized.bufferSize());
            quadraturePoints.front().evaluate(
                source,
                sourceResult,
                sourceBuffer);
            quadraturePoints.front().evaluate(
                deserialized,
                deserializedResult,
                deserializedBuffer);
            CIE_TEST_REQUIRE(std::any_of(
                sourceResult.begin(),
                sourceResult.end(),
                [] (Scalar value) { return std::abs(value) > 1e-12; }));
            for (std::size_t iValue=0ul; iValue<sourceResult.size(); ++iValue)
                CIE_TEST_CHECK(deserializedResult[iValue] == Approx(sourceResult[iValue]).margin(1e-12));
        }

        const auto makeSink = [] (std::vector<std::vector<Scalar>>* pBatches) {
            return [pBatches] (
                std::span<const VertexID> cellIDs,
                std::span<const Scalar> values) {
                    CIE_TEST_CHECK(cellIDs.size() == 2ul);
                    CIE_TEST_CHECK(values.size() == 2ul * 16ul);
                    pBatches->emplace_back(values.begin(), values.end());
            };
        };

        const IntegrandProcessor<Dimension,Integrand,Scalar>::Properties properties {
            .integrandBatchSize = 8ul,
            .verbosity = 0u};

        // Reference results from the host processor.
        std::vector<std::vector<Scalar>> hostBatches;
        {
            IntegrandProcessor<Dimension,Integrand,Scalar> hostProcessor;
            hostProcessor.process(
                std::span<const ActiveTestCell2D>(cells),
                quadratureRuleFactory,
                integrandFactory,
                makeSink(&hostBatches),
                properties);
        }
        CIE_TEST_REQUIRE(hostBatches.size() == 1ul);
        CIE_TEST_REQUIRE(hostBatches.front().size() == 2ul * 16ul);
        CIE_TEST_REQUIRE(std::any_of(
            hostBatches.front().begin(),
            hostBatches.front().end(),
            [] (Scalar value) { return std::abs(value) > 1e-12; }));

        // Results from the SYCL processor, which serializes the integrands, moves
        // them to the device, and deserializes them into a device memory pool.
        std::vector<std::vector<Scalar>> syclBatches;
        {
            sycl::queue queue(sycl::default_selector_v);
            SYCLIntegrandProcessor<Dimension,Integrand,Scalar> syclProcessor(
                std::make_shared<sycl::queue>(std::move(queue)),
                &stats);
            syclProcessor.process(
                std::span<const ActiveTestCell2D>(cells),
                quadratureRuleFactory,
                integrandFactory,
                makeSink(&syclBatches),
                properties);
        }
        CIE_TEST_REQUIRE(syclBatches.size() == 1ul);
        CIE_TEST_REQUIRE(syclBatches.front().size() == hostBatches.front().size());
        for (std::size_t iComponent=0ul; iComponent<syclBatches.front().size(); ++iComponent) {
            INFO("component " << iComponent << ", SYCL=" << syclBatches.front()[iComponent]
                << ", host=" << hostBatches.front()[iComponent]);
            CIE_TEST_CHECK(syclBatches.front()[iComponent] == Approx(hostBatches.front()[iComponent]).margin(1e-8));
        }
    }
#endif // CIE_ENABLE_SYCL


} // namespace cie::fem
