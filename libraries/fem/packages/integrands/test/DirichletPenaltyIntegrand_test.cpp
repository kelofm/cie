// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


struct DirichletPenaltyTest : maths::ExpressionTraits<double> {
    using ExpressionTraits<double>::ConstSpan;
    using ExpressionTraits<double>::Span;
    using ExpressionTraits<double>::BufferSpan;

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
        spatialTransform);
    CIE_TEST_REQUIRE(integrand.size() == 4 * 4 + 4);
    CIE_TEST_CHECK(integrand.bufferSize() == 4 + 2 + 2 + 1 + 5);

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
    } // for rSamplePoint, rReferences in references
}


} // namespace cie::fem
