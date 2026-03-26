// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"


namespace cie::fem {


CIE_TEST_CASE("LinearIsotropicStiffnessIntegrand", "[integrands]") {
    CIE_TEST_CASE_INIT("LinearIsotropicStiffnessIntegrand")
    using Scalar = double;
    constexpr unsigned Dimension = 2u;

    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;

    // Define a bilinear ansatz space.
    const auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
    });

    // Compute the derivatives of the ansatz space.
    const auto pAnsatzDerivatives = std::make_shared<Ansatz::Derivative>(pAnsatzSpace->makeDerivative());

    // Construct the integrand without a buffer.
    constexpr Scalar modulus = 10.0;
    LinearIsotropicStiffnessIntegrand<Ansatz::Derivative> integrand(modulus, *pAnsatzDerivatives);
    CIE_TEST_REQUIRE(integrand.size() == 16);
    CIE_TEST_REQUIRE(integrand.bufferSize() == 8);

    // Set buffer.
    std::vector<Scalar> stiffness(integrand.size());
    std::vector<Scalar> buffer(integrand.bufferSize());

    // Check stiffness values.
    std::vector<Scalar> reference(integrand.size());
    StaticArray<Scalar,Dimension> location;

    DynamicArray<std::pair<
        StaticArray<Scalar,Dimension>,  //< sample point
        std::vector<Scalar>             //< reference values
    >> references;
    references.reserve(4);

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, -1.0}},
        std::vector<Scalar> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.25,  0.00, -0.25,
             0.00,  0.00,  0.25, -0.25,
             0.00, -0.25, -0.25,  0.50
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0, -1.0}},
        std::vector<Scalar> {{
             0.25,  0.00, -0.25,  0.00,
             0.00,  0.00,  0.00,  0.00,
            -0.25,  0.00,  0.50, -0.25,
             0.00,  0.00, -0.25,  0.25
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, 1.0}},
        std::vector<Scalar> {{
             0.25, -0.25,  0.00,  0.00,
            -0.25,  0.50,  0.00, -0.25,
             0.00,  0.00,  0.00,  0.00,
             0.00, -0.25,  0.00,  0.25
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0,  1.0}},
        std::vector<Scalar> {{
             0.50, -0.25, -0.25,  0.00,
            -0.25,  0.25,  0.00,  0.00,
            -0.25,  0.00,  0.25,  0.00,
             0.00,  0.00,  0.00,  0.00
        }}
    );

    for (const auto& [rSamplePoint, rReference] : references) {
        std::vector<Scalar> result(integrand.size());
        CIE_TEST_CHECK_NOTHROW(integrand.evaluate(rSamplePoint, result, buffer));
        for (unsigned iComponent=0u; iComponent<rReference.size(); ++iComponent) {
            CIE_TEST_CHECK(result[iComponent] == Approx(modulus * rReference[iComponent]).margin(1e-14));
        } // for iComponent in range(rReference.size())
    } // for rSamplePoint, rReferences in references
}


} // namespace cie::fem
