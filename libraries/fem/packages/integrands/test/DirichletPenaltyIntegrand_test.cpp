// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


struct DirichletPenaltyTest
    :   maths::ExpressionTraits<double>,
        cie::io::TriviallySerializableBase {
    using ExpressionTraits<double>::ConstSpan;
    using ExpressionTraits<double>::Span;
    using ExpressionTraits<double>::BufferSpan;

    template <template <class ...> class, class>
    using Rebind = DirichletPenaltyTest;

    void evaluate(ConstSpan in, Span out, BufferSpan buffer) const {
        CIE_TEST_CHECK(in.size() == 2);
        CIE_TEST_CHECK(out.size() == this->size());
        CIE_TEST_CHECK(buffer.size() == this->bufferSize());
        out[0] = in[0] + 2 * in[1];
    }

    unsigned size() const noexcept {
        return 1u;
    }

    unsigned bufferSize() const noexcept {
        return 5u;
    }
}; // DirichletPenaltyTest


CIE_TEST_CASE("DirichletPenaltyIntegrand", "[integrands]")
{
    CIE_TEST_CASE_INIT("DirichletPenaltyIntegrand")
    using Scalar = double;
    constexpr unsigned Dimension = 2u;

    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;

    // Define a bilinear ansatz space.
    const auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
    });

    // Construct the integrand without a buffer.
    constexpr Scalar penalty = 10.0;
    const maths::IdentityTransform<Scalar,Dimension> spatialTransform;
    const DirichletPenaltyTest dirichlet;
    DirichletPenaltyIntegrand<
        DirichletPenaltyTest,
        Ansatz,
        maths::IdentityTransform<Scalar,Dimension>,
        CellBase<Dimension,Scalar,maths::IdentityTransform<Scalar,Dimension>>> integrand(
            dirichlet,
            penalty,
            *pAnsatzSpace,
            spatialTransform,
            spatialTransform.makeInverse(),
            spatialTransform.makeDerivative());
    CIE_TEST_REQUIRE(integrand.size() == 4 * 4 + 4);
    CIE_TEST_CHECK(integrand.bufferSize() == 4 + 2 + 2 + 1 + 5);

    using Cell = CellBase<Dimension,Scalar,maths::IdentityTransform<Scalar,Dimension>>;
    DirichletPenaltyLHSIntegrand<Ansatz,maths::IdentityTransform<Scalar,Dimension>,Cell> lhsIntegrand(
        penalty,
        *pAnsatzSpace,
        spatialTransform,
        spatialTransform.makeInverse(),
        spatialTransform.makeDerivative());
    DirichletPenaltyRHSIntegrand<
        DirichletPenaltyTest,
        Ansatz,
        maths::IdentityTransform<Scalar,Dimension>,
        Cell> rhsIntegrand(
            dirichlet,
            penalty,
            *pAnsatzSpace,
            spatialTransform,
            spatialTransform.makeInverse(),
            spatialTransform.makeDerivative());
    static_assert(maths::SliceableExpression<decltype(lhsIntegrand)>);
    static_assert(maths::SliceableExpression<decltype(rhsIntegrand)>);
    CIE_TEST_REQUIRE(lhsIntegrand.size() == 16);
    CIE_TEST_REQUIRE(rhsIntegrand.size() == 4);

    // Set buffers.
    std::vector<Scalar> output;
    output.resize(integrand.size());
    std::vector<Scalar> buffer;
    buffer.resize(integrand.bufferSize());

    // Check mass values.
    DynamicArray<std::pair<
        StaticArray<Scalar,Dimension>,  //< sample point
        StaticArray<Scalar,20>          //< reference values
    >> references;
    references.reserve(4);

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, -1.0}},
        StaticArray<Scalar,20> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  1.00,

             0.00,  0.00,  0.00, -3.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0, -1.0}},
        StaticArray<Scalar,20> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  1.00,  0.00,
             0.00,  0.00,  0.00,  0.00,

             0.00,  0.00, -1.00,  0.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, 1.0}},
        StaticArray<Scalar,20> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  1.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,

             0.00,  1.00,  0.00,  0.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0,  1.0}},
        StaticArray<Scalar,20> {{
             1.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,

             3.00,  0.00,  0.00,  0.00
        }}
    );

    for (const auto& [rSamplePoint, rReference] : references) {
        std::vector<Scalar> result;
        result.resize(integrand.size());
        CIE_TEST_CHECK_NOTHROW(integrand.evaluate(rSamplePoint, result, buffer));
        for (unsigned iComponent=0u; iComponent<rReference.size(); ++iComponent) {
            CIE_TEST_CHECK(result[iComponent] == Approx(penalty * rReference[iComponent]).margin(1e-14));
        } // for iComponent in range(rReference.size())

        std::vector<Scalar> lhsResult(lhsIntegrand.size());
        std::vector<Scalar> lhsBuffer(lhsIntegrand.bufferSize());
        CIE_TEST_CHECK_NOTHROW(lhsIntegrand.evaluate(rSamplePoint, lhsResult, lhsBuffer));
        for (std::size_t iComponent=0ul; iComponent<lhsResult.size(); ++iComponent)
            CIE_TEST_CHECK(lhsResult[iComponent] == Approx(result[iComponent]).margin(1e-14));

        std::vector<Scalar> rhsResult(rhsIntegrand.size());
        std::vector<Scalar> rhsBuffer(rhsIntegrand.bufferSize());
        CIE_TEST_CHECK_NOTHROW(rhsIntegrand.evaluate(rSamplePoint, rhsResult, rhsBuffer));
        for (std::size_t iComponent=0ul; iComponent<rhsResult.size(); ++iComponent)
            CIE_TEST_CHECK(rhsResult[iComponent] == Approx(result[lhsResult.size() + iComponent]).margin(1e-14));

        // Assemble the LHS from row/column tile pairs.
        std::vector<Ansatz::Context> tiling;
        tiling.emplace_back(std::array<Ansatz::Context::Value,2 * Dimension> {0, 1, 0, 2});
        tiling.emplace_back(std::array<Ansatz::Context::Value,2 * Dimension> {1, 2, 0, 2});
        const std::array<Ansatz::Context::Value,4> lhsOutputBounds {2, 2, 2, 2};
        std::vector<Scalar> tiledLHS(lhsIntegrand.size(), Scalar(0));
        for (const auto& rRowTile : tiling) {
            for (const auto& rColumnTile : tiling) {
                const decltype(lhsIntegrand)::InContext inContext{rRowTile,rColumnTile};
                std::vector<Scalar> tileResult(lhsIntegrand.size(inContext));
                std::vector<Scalar> tileBuffer(lhsIntegrand.bufferSize(inContext));
                CIE_TEST_CHECK_NOTHROW(lhsIntegrand.evaluate(rSamplePoint, tileResult, tileBuffer, inContext));
                const auto outContext = lhsIntegrand.outContext(inContext);
                std::vector<Ansatz::Context::Value> outputIndices(tileResult.size());
                CIE_TEST_REQUIRE(outContext.count() == tileResult.size());
                CIE_TEST_REQUIRE_NOTHROW(outContext.expand(lhsOutputBounds, outputIndices));
                for (std::size_t iComponent=0ul; iComponent<tileResult.size(); ++iComponent)
                    tiledLHS[outputIndices[iComponent]] += tileResult[iComponent];
            }
        }
        for (std::size_t iComponent=0ul; iComponent<tiledLHS.size(); ++iComponent)
            CIE_TEST_CHECK(tiledLHS[iComponent] == Approx(lhsResult[iComponent]).margin(1e-14));

        // Assemble the RHS from single-context tiles.
        const std::array<Ansatz::Context::Value,3> rhsOutputBounds {2, 2, 1};
        std::vector<Scalar> tiledRHS(rhsIntegrand.size(), Scalar(0));
        for (const auto& rTile : tiling) {
            const decltype(rhsIntegrand)::InContext inContext{rTile};
            std::vector<Scalar> tileResult(rhsIntegrand.size(inContext));
            std::vector<Scalar> tileBuffer(rhsIntegrand.bufferSize(inContext));
            CIE_TEST_CHECK_NOTHROW(rhsIntegrand.evaluate(rSamplePoint, tileResult, tileBuffer, inContext));
            const auto outContext = rhsIntegrand.outContext(inContext);
            std::vector<Ansatz::Context::Value> outputIndices(tileResult.size());
            CIE_TEST_REQUIRE(outContext.count() == tileResult.size());
            CIE_TEST_REQUIRE_NOTHROW(outContext.expand(rhsOutputBounds, outputIndices));
            for (std::size_t iComponent=0ul; iComponent<tileResult.size(); ++iComponent)
                tiledRHS[outputIndices[iComponent]] += tileResult[iComponent];
        }
        for (std::size_t iComponent=0ul; iComponent<tiledRHS.size(); ++iComponent)
            CIE_TEST_CHECK(tiledRHS[iComponent] == Approx(rhsResult[iComponent]).margin(1e-14));
    } // for rSamplePoint, rReferences in references
}


} // namespace cie::fem
