// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/IndentityTransform.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


struct DirichletPenaltyTest : maths::ExpressionTraits<double>
{
    using ExpressionTraits<double>::ConstSpan;
    using ExpressionTraits<double>::Span;

    void evaluate(ConstSpan in, Span out) const {
        CIE_TEST_CHECK(in.size() == 2);
        CIE_TEST_CHECK(out.size() == 1);
        out[0] = in[0] + 2 * in[1];
    }

    unsigned size() const noexcept {
        return 1u;
    }
}; // DirichletPenaltyTest


CIE_TEST_CASE("DirichletPenaltyIntegrand", "[integrands]")
{
    CIE_TEST_CASE_INIT("DirichletPenaltyIntegrand")
    using Scalar = double;
    constexpr unsigned Dimension = 2u;

    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;
    mp::ThreadPoolBase threadPool;

    // Define a bilinear ansatz space.
    const auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
    }, threadPool);

    // Construct the integrand without a buffer.
    constexpr Scalar penalty = 10.0;
    const maths::IdentityTransform<Scalar,Dimension> spatialTransform;
    const DirichletPenaltyTest dirichlet;
    DirichletPenaltyIntegrand<DirichletPenaltyTest,Ansatz,maths::IdentityTransform<Scalar,Dimension>> integrand(
        dirichlet,
        penalty,
        *pAnsatzSpace,
        spatialTransform);
    CIE_TEST_CHECK(integrand.size() == 4 * 4 + 4);
    CIE_TEST_CHECK(integrand.getMinBufferSize() == 4 + 2 + 1);

    // Set buffer.
    StaticArray<Scalar,7> buffer;
    StaticArray<Scalar,20> output;

    #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        StaticArray<Scalar,Dimension> dummy;
        std::fill(dummy.begin(), dummy.end(), 0);

        // Attempt to evaluate without setting a buffer.
        CIE_TEST_CHECK_THROWS(integrand.evaluate(dummy, {output.data(), output.data() + output.size()}));

        // Attempt to set insufficiently sized buffers.
        CIE_TEST_CHECK_THROWS(integrand.setBuffer({}));

        for (unsigned bufferSize : {0u, 1u, 6u}) {
            CIE_TEST_CHECK_THROWS(integrand.setBuffer({buffer.data(), bufferSize}));
        }

        // Attempt to construct with insufficiently sized buffers.
        for (unsigned bufferSize : {0u, 1u, 6u}) {
            CIE_TEST_CHECK_THROWS(integrand = DirichletPenaltyIntegrand<Ansatz::Derivative>(
                1.0,
                pAnsatzSpace,
                {buffer.data(), bufferSize}));
        }
    #endif

    CIE_TEST_CHECK_NOTHROW(integrand.setBuffer(buffer));

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
        StaticArray<Scalar,20> result;
        CIE_TEST_CHECK_NOTHROW(integrand.evaluate(rSamplePoint, result));
        for (unsigned iComponent=0u; iComponent<rReference.size(); ++iComponent) {
            CIE_TEST_CHECK(result[iComponent] == Approx(penalty * rReference[iComponent]).margin(1e-14));
        } // for iComponent in range(rReference.size())
    } // for rSamplePoint, rReferences in references
}


} // namespace cie::fem
