// --- Internal Includes ---
#include "packages/maths/inc/LegendrePolynomial.hpp"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <ranges>


namespace cie::fem::maths {


CIE_TEST_CASE("LegendrePolynomial", "[maths]")
{
    CIE_TEST_CASE_INIT("LegendrePolynomial")

    using P = LegendrePolynomial<double>;

    //for (unsigned p=0u; p<10; ++p) {
    //    CIE_TEST_REQUIRE_NOTHROW(P(p));
    //    P polynomial(p);
    //    std::cout << p << ": ";
    //    std::ranges::copy(
    //        std::views::transform(
    //            polynomial.coefficients(),
    //            [] (double c) {return std::to_string(c);}),
    //        std::ostream_iterator<std::string>(std::cout, " ")
    //    );
    //    std::cout << "\n";
    //}

    {
        CIE_TEST_CASE_INIT("polynomial order 0")
        CIE_TEST_REQUIRE_NOTHROW(P(0));
        P p(0);
        CIE_TEST_REQUIRE(p.coefficients().size() == 1);
        CIE_TEST_CHECK(p.coefficients()[0] == Approx(1.0).margin(1e-12));
    }

    {
        CIE_TEST_CASE_INIT("polynomial order 1")
        CIE_TEST_REQUIRE_NOTHROW(P(1));
        P p(1);
        CIE_TEST_REQUIRE(p.coefficients().size() == 2);
        CIE_TEST_CHECK(p.coefficients()[0] == Approx(0.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[1] == Approx(1.0).margin(1e-12));
    }

    {
        CIE_TEST_CASE_INIT("polynomial order 2")
        CIE_TEST_REQUIRE_NOTHROW(P(2));
        P p(2);
        CIE_TEST_REQUIRE(p.coefficients().size() == 3);
        CIE_TEST_CHECK(p.coefficients()[0] == Approx(-0.5).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[1] == Approx( 0.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[2] == Approx( 1.5).margin(1e-12));
    }

    {
        CIE_TEST_CASE_INIT("polynomial order 3")
        CIE_TEST_REQUIRE_NOTHROW(P(3));
        P p(3);
        CIE_TEST_REQUIRE(p.coefficients().size() == 4);
        CIE_TEST_CHECK(p.coefficients()[0] == Approx(-0.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[1] == Approx(-1.5).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[2] == Approx( 0.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[3] == Approx( 2.5).margin(1e-12));
    }

    {
        CIE_TEST_CASE_INIT("polynomial order 10")
        CIE_TEST_REQUIRE_NOTHROW(P(10));
        P p(10);
        CIE_TEST_REQUIRE(p.coefficients().size() == 11);
        CIE_TEST_CHECK(p.coefficients()[ 0] == Approx(-    63.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 1] == Approx(      0.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 2] == Approx(   3465.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 3] == Approx(      0.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 4] == Approx(- 30030.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 5] == Approx(      0.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 6] == Approx(  90090.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 7] == Approx(      0.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 8] == Approx(-109395.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[ 9] == Approx(      0.0 / 256.0).margin(1e-12));
        CIE_TEST_CHECK(p.coefficients()[10] == Approx(  46189.0 / 256.0).margin(1e-12));
    }
}


} // namespace cie::fem::maths
