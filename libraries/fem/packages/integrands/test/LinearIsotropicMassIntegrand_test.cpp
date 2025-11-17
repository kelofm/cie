// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/integrands/inc/LinearIsotropicMassIntegrand.hpp"


namespace cie::fem {


CIE_TEST_CASE("LinearIsotropicMassIntegrand", "[integrands]")
{
    CIE_TEST_CASE_INIT("LinearIsotropicMassIntegrand")
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
    constexpr Scalar modulus = 10.0;
    LinearIsotropicMassIntegrand<Ansatz> integrand(modulus, *pAnsatzSpace);
    CIE_TEST_CHECK(integrand.size() == 16);
    CIE_TEST_CHECK(integrand.getMinBufferSize() == 4);

    // Set buffer.
    StaticArray<Scalar,4> buffer;
    StaticArray<Scalar,16> mass;

    #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
        StaticArray<Scalar,Dimension> dummy;
        std::fill(dummy.begin(), dummy.end(), 0);

        // Attempt to evaluate without setting a buffer.
        CIE_TEST_CHECK_THROWS(integrand.evaluate(dummy, {mass.data(), mass.data() + mass.size()}));

        // Attempt to set insufficiently sized buffers.
        CIE_TEST_CHECK_THROWS(integrand.setBuffer({}));

        for (unsigned bufferSize : {0u, 1u, 3u}) {
            CIE_TEST_CHECK_THROWS(integrand.setBuffer({buffer.data(), bufferSize}));
        }

        // Attempt to construct with insufficiently sized buffers.
        for (unsigned bufferSize : {0u, 1u, 3u}) {
            CIE_TEST_CHECK_THROWS(integrand = LinearIsotropicMassIntegrand<Ansatz::Derivative>(
                1.0,
                pAnsatzSpace,
                {buffer.data(), bufferSize}));
        }
    #endif

    CIE_TEST_CHECK_NOTHROW(integrand.setBuffer({buffer.data(), buffer.size()}));

    // Check mass values.
    DynamicArray<std::pair<
        StaticArray<Scalar,Dimension>,  //< sample point
        StaticArray<Scalar,16>          //< reference values
    >> references;
    references.reserve(4);

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, -1.0}},
        StaticArray<Scalar,16> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  1.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0, -1.0}},
        StaticArray<Scalar,16> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  1.00,  0.00,
             0.00,  0.00,  0.00,  0.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{-1.0, 1.0}},
        StaticArray<Scalar,16> {{
             0.00,  0.00,  0.00,  0.00,
             0.00,  1.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00
        }}
    );

    references.emplace_back(
        StaticArray<Scalar,Dimension> {{ 1.0,  1.0}},
        StaticArray<Scalar,16> {{
             1.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00,
             0.00,  0.00,  0.00,  0.00
        }}
    );

    for (const auto& [rSamplePoint, rReference] : references) {
        StaticArray<Scalar,16> result;
        CIE_TEST_CHECK_NOTHROW(integrand.evaluate(rSamplePoint, result));
        for (unsigned iComponent=0u; iComponent<rReference.size(); ++iComponent) {
            CIE_TEST_CHECK(result[iComponent] == Approx(modulus * rReference[iComponent]).margin(1e-14));
        } // for iComponent in range(rReference.size())
    } // for rSamplePoint, rReferences in references
}


} // namespace cie::fem
