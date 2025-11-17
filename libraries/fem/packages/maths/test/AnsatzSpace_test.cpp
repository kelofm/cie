// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/kernel.hpp"

namespace cie::fem::maths {


CIE_TEST_CASE( "AnsatzSpace", "[maths]" )
{
    CIE_TEST_CASE_INIT( "AnsatzSpace" )

    using Basis = Polynomial<double>;
    using AnsatzSpace = AnsatzSpace<Basis,2>;

    DynamicArray<Basis> basisFunctions {
        Basis(Basis::Coefficients { 0.5, -0.5}),
        Basis(Basis::Coefficients { 0.5,  0.5}),
        Basis(Basis::Coefficients { 0.0,  0.0, 1.0})
    };

    CIE_TEST_REQUIRE_NOTHROW(AnsatzSpace(basisFunctions));
    AnsatzSpace ansatzSpace(basisFunctions);
    CIE_TEST_REQUIRE(ansatzSpace.size() == 9);

    using Point = Kernel<2,double>::Point;
    double x = 11.0;
    double y = 9.0;
    Point p0 {-x, -y};
    Point p1 {-x,  y};
    Point p2 { x, -y};
    Point p3 { x,  y};

    StaticArray<double,9> results;

    CIE_TEST_CHECK_NOTHROW(ansatzSpace.evaluate(p0, results));
    CIE_TEST_CHECK(results[0] == Approx(  30.0));
    CIE_TEST_CHECK(results[1] == Approx( -25.0));
    CIE_TEST_CHECK(results[2] == Approx( 605.0));
    CIE_TEST_CHECK(results[3] == Approx( -24.0));
    CIE_TEST_CHECK(results[4] == Approx(  20.0));
    CIE_TEST_CHECK(results[5] == Approx(-484.0));
    CIE_TEST_CHECK(results[6] == Approx( 486.0));
    CIE_TEST_CHECK(results[7] == Approx(-405.0));
    CIE_TEST_CHECK(results[8] == Approx(9801.0));

    CIE_TEST_CHECK_NOTHROW(ansatzSpace.evaluate(p1, results));
    CIE_TEST_CHECK(results[0] == Approx( -24.0));
    CIE_TEST_CHECK(results[1] == Approx(  20.0));
    CIE_TEST_CHECK(results[2] == Approx(-484.0));
    CIE_TEST_CHECK(results[3] == Approx(  30.0));
    CIE_TEST_CHECK(results[4] == Approx( -25.0));
    CIE_TEST_CHECK(results[5] == Approx( 605.0));
    CIE_TEST_CHECK(results[6] == Approx( 486.0));
    CIE_TEST_CHECK(results[7] == Approx(-405.0));
    CIE_TEST_CHECK(results[8] == Approx(9801.0));

    CIE_TEST_CHECK_NOTHROW(ansatzSpace.evaluate(p2, results));
    CIE_TEST_CHECK(results[0] == Approx( -25.0));
    CIE_TEST_CHECK(results[1] == Approx(  30.0));
    CIE_TEST_CHECK(results[2] == Approx( 605.0));
    CIE_TEST_CHECK(results[3] == Approx(  20.0));
    CIE_TEST_CHECK(results[4] == Approx( -24.0));
    CIE_TEST_CHECK(results[5] == Approx(-484.0));
    CIE_TEST_CHECK(results[6] == Approx(-405.0));
    CIE_TEST_CHECK(results[7] == Approx( 486.0));
    CIE_TEST_CHECK(results[8] == Approx(9801.0));

    CIE_TEST_CHECK_NOTHROW(ansatzSpace.evaluate(p3, results));
    CIE_TEST_CHECK(results[0] == Approx(  20.0));
    CIE_TEST_CHECK(results[1] == Approx( -24.0));
    CIE_TEST_CHECK(results[2] == Approx(-484.0));
    CIE_TEST_CHECK(results[3] == Approx( -25.0));
    CIE_TEST_CHECK(results[4] == Approx(  30.0));
    CIE_TEST_CHECK(results[5] == Approx( 605.0));
    CIE_TEST_CHECK(results[6] == Approx(-405.0));
    CIE_TEST_CHECK(results[7] == Approx( 486.0));
    CIE_TEST_CHECK(results[8] == Approx(9801.0));
}


CIE_TEST_CASE( "AnsatzSpaceDerivative", "[maths]" )
{
    CIE_TEST_CASE_INIT( "AnsatzSpaceDerivative" )

    using Basis = Polynomial<double>;
    using AnsatzSpace = AnsatzSpace<Basis,2>;

    DynamicArray<Basis> basisFunctions {
        Basis(Basis::Coefficients {0.5, -0.5}),
        Basis(Basis::Coefficients {0.5, 0.5})
    };

    CIE_TEST_REQUIRE_NOTHROW(AnsatzSpace(basisFunctions).makeDerivative());
    const auto ansatzDerivative = AnsatzSpace(basisFunctions).makeDerivative();
    CIE_TEST_REQUIRE(ansatzDerivative.size() == 8);

    using Point = Kernel<2,double>::Point;
    const double x = 2;
    const double y = 3;
    const Point p0 {-x, -y};
    const Point p1 { x, -y};
    const Point p2 {-x,  y};
    const Point p3 { x,  y};

    Eigen::Matrix<double,4,2> buffer;

    ansatzDerivative.evaluate(p0, {buffer.data(), buffer.data() + buffer.size()});
    CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));
    CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));
    CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));
    CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));
    CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));
    CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));
    CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));
    CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));

    ansatzDerivative.evaluate(p1, {buffer.data(), buffer.data() + buffer.size()});
    CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));
    CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));
    CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));
    CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));
    CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));
    CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));
    CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));
    CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));

    ansatzDerivative.evaluate(p2, {buffer.data(), buffer.data() + buffer.size()});
    CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));
    CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));
    CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));
    CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));
    CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));
    CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));
    CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));
    CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));

    ansatzDerivative.evaluate(p3, {buffer.data(), buffer.data() + buffer.size()});
    CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));
    CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));
    CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));
    CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));
    CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));
    CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));
    CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));
    CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));
}


} // namespace cie::fem::maths
