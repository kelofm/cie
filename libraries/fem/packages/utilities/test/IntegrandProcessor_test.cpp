// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/utilities/inc/IntegrandProcessor.hpp"


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


} // namespace cie::fem
