// ---- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"

// --- STL Includes ---
#include <cmath>


namespace cie::fem {


CIE_TEST_CASE( "GaussLegendreQuadrature", "[numeric]" )
{
    CIE_TEST_CASE_INIT( "GaussLegendreQuadrature" )

    const double maxAbsoluteError = 5e-15;
    const Size maxIterations = 20;

    {
        CIE_TEST_CASE_INIT("order = 5")

        const Size integrationOrder = 5;
        std::pair<std::vector<double>,std::vector<double>> reference {
            {
                -0.9061798459386640,
                -0.5384693101056831,
                0.0000000000000000,
                0.5384693101056831,
                0.9061798459386640
            },
            {
                0.2369268850561894,
                0.4786286704993662,
                0.5688888888888890,
                0.4786286704993662,
                0.2369268850561894
            }
        };

        GaussLegendreQuadrature<double> quadrature(integrationOrder, maxAbsoluteError, maxIterations);
        CIE_TEST_REQUIRE(quadrature.nodes().size() == integrationOrder);
        CIE_TEST_REQUIRE(quadrature.weights().size() == integrationOrder);

        for (Size index=0; index<integrationOrder; ++index) {
            CIE_TEST_CHECK(quadrature.nodes()[index] == Approx(reference.first[index]).margin(maxAbsoluteError));
            CIE_TEST_CHECK(quadrature.weights()[index] == Approx(reference.second[index]).margin(maxAbsoluteError));
        }
    }

    {
        CIE_TEST_CASE_INIT("order = 10")

        const Size integrationOrder = 10;
        std::pair<std::vector<double>,std::vector<double>> reference {
            {
                -0.9739065285171717,
                -0.8650633666889845,
                -0.6794095682990244,
                -0.4333953941292472,
                -0.1488743389816312,
                0.1488743389816312,
                0.4333953941292472,
                0.6794095682990244,
                0.8650633666889845,
                0.9739065285171717
            },
            {
                0.0666713443086881,
                0.1494513491505804,
                0.2190863625159820,
                0.2692667193099965,
                0.2955242247147530,
                0.2955242247147530,
                0.2692667193099965,
                0.2190863625159820,
                0.1494513491505804,
                0.0666713443086881,
            }
        };

        GaussLegendreQuadrature<double> quadrature(integrationOrder, maxAbsoluteError, maxIterations);
        CIE_TEST_REQUIRE(quadrature.nodes().size() == integrationOrder);
        CIE_TEST_REQUIRE(quadrature.weights().size() == integrationOrder);

        for (Size index=0; index<integrationOrder; ++index)
        {
            CIE_TEST_CHECK(quadrature.nodes()[index] == Approx(reference.first[index]).margin(maxAbsoluteError));
            CIE_TEST_CHECK(quadrature.weights()[index] == Approx(reference.second[index]).margin(maxAbsoluteError));
        }
    }

    {
        CIE_TEST_CASE_INIT("order = 20")

        const Size integrationOrder = 20;
        std::pair<std::vector<double>,std::vector<double>> reference {
            {
                -0.9931285991850949,
                -0.9639719272779138,
                -0.9122344282513258,
                -0.8391169718222188,
                -0.7463319064601508,
                -0.6360536807265150,
                -0.5108670019508271,
                -0.3737060887154195,
                -0.2277858511416451,
                -0.0765265211334973,
                0.0765265211334973,
                0.2277858511416451,
                0.3737060887154195,
                0.5108670019508271,
                0.6360536807265150,
                0.7463319064601508,
                0.8391169718222188,
                0.9122344282513258,
                0.9639719272779138,
                0.9931285991850949
            },
            {
                0.0176140071391533,
                0.0406014298003862,
                0.0626720483341094,
                0.0832767415767047,
                0.1019301198172403,
                0.1181945319615182,
                0.1316886384491765,
                0.1420961093183819,
                0.1491729864726037,
                0.1527533871307258,
                0.1527533871307258,
                0.1491729864726037,
                0.1420961093183819,
                0.1316886384491765,
                0.1181945319615182,
                0.1019301198172403,
                0.0832767415767047,
                0.0626720483341094,
                0.0406014298003862,
                0.0176140071391533
            }
        };

        GaussLegendreQuadrature<double> quadrature(integrationOrder, maxAbsoluteError, maxIterations);
        CIE_TEST_REQUIRE(quadrature.nodes().size() == integrationOrder);
        CIE_TEST_REQUIRE(quadrature.weights().size() == integrationOrder);

        for (Size index=0; index<integrationOrder; ++index) {
            CIE_TEST_CHECK(quadrature.nodes()[index] == Approx(reference.first[index]).margin(maxAbsoluteError));
            CIE_TEST_CHECK(quadrature.weights()[index] == Approx(reference.second[index]).margin(maxAbsoluteError));
        }
    }
}


} // namespace cie::fem