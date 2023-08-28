#include "packages/matrix/inc/StaticEigenMatrix.hpp"
#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- LinAlg Includes ---
#include "packages/spatial/inc/Quaternion.hpp"
#include "packages/overloads/inc/matrix_operators.hpp"

// --- STL Includes ---
#include <memory>

namespace cie::linalg {



CIE_TEST_CASE("Quaternion", "[spatial]")
{
    std::unique_ptr<Quaternion<double>> p_quaternion;

    CIE_TEST_REQUIRE_NOTHROW(p_quaternion = std::make_unique<Quaternion<double>>());
    CIE_TEST_REQUIRE(p_quaternion->components().size() == 4);
    CIE_TEST_CHECK((*p_quaternion)[0] == Approx(1.0));
    CIE_TEST_CHECK((*p_quaternion)[1] == Approx(0.0));
    CIE_TEST_CHECK((*p_quaternion)[2] == Approx(0.0));
    CIE_TEST_CHECK((*p_quaternion)[3] == Approx(0.0));

    CIE_TEST_REQUIRE_NOTHROW(p_quaternion.reset(new Quaternion<double>({0.0,1.0,2.0,3.0})));
    CIE_TEST_REQUIRE(p_quaternion->components().size() == 4);
    CIE_TEST_CHECK((*p_quaternion)[0] == Approx(0.0));
    CIE_TEST_CHECK((*p_quaternion)[1] == Approx(1.0));
    CIE_TEST_CHECK((*p_quaternion)[2] == Approx(2.0));
    CIE_TEST_CHECK((*p_quaternion)[3] == Approx(3.0));

    CIE_TEST_REQUIRE_NOTHROW(p_quaternion.reset(new Quaternion<double>(*p_quaternion)));
    CIE_TEST_REQUIRE(p_quaternion->components().size() == 4);
    CIE_TEST_CHECK((*p_quaternion)[0] == Approx(0.0));
    CIE_TEST_CHECK((*p_quaternion)[1] == Approx(1.0));
    CIE_TEST_CHECK((*p_quaternion)[2] == Approx(2.0));
    CIE_TEST_CHECK((*p_quaternion)[3] == Approx(3.0));

    Quaternion<double> p = *p_quaternion;
    CIE_TEST_REQUIRE(p.components().size() == 4);
    CIE_TEST_CHECK(p[0] == Approx(0.0));
    CIE_TEST_CHECK(p[1] == Approx(1.0));
    CIE_TEST_CHECK(p[2] == Approx(2.0));
    CIE_TEST_CHECK(p[3] == Approx(3.0));

    CIE_TEST_CHECK(p.begin() == p.components().begin());
    CIE_TEST_CHECK(p.end() == p.components().end());

    CIE_TEST_CHECK(p.normSquared() == Approx(14.0));
    CIE_TEST_CHECK_NOTHROW(p.normalize());
    CIE_TEST_CHECK(p.normSquared() == Approx(1.0));
    const double norm = std::sqrt(14.0);
    CIE_TEST_CHECK(p[0] == Approx(0.0));
    CIE_TEST_CHECK(p[1] == Approx(1.0 / norm));
    CIE_TEST_CHECK(p[2] == Approx(2.0 / norm));
    CIE_TEST_CHECK(p[3] == Approx(3.0 / norm));

    // Check construction from axis + angle
    StaticArray<double,3> axis { 0.0, 0.0, 1.0 };
    const double angle = 90.0 * M_PI / 180.0;
    CIE_TEST_CHECK_NOTHROW(p.loadFromAxisAndAngle(axis, angle));
    const double sqrt2Over2 = std::sqrt(2.0) / 2.0;
    CIE_TEST_CHECK(p[0] == Approx(sqrt2Over2));
    CIE_TEST_CHECK(p[1] == Approx(0.0));
    CIE_TEST_CHECK(p[2] == Approx(0.0));
    CIE_TEST_CHECK(p[3] == Approx(sqrt2Over2));

    StaticEigenMatrix<double,3,3> rotationMatrix;
    CIE_TEST_CHECK_NOTHROW(p.toRotationMatrix(rotationMatrix));

    StaticEigenMatrix<double,3,1> vector;
    vector(0, 0) = 1.0;
    vector(1, 0) = 0.0;
    vector(2, 0) = 0.0;
    CIE_TEST_CHECK_NOTHROW(vector = rotationMatrix * vector);
    CIE_TEST_CHECK(vector(0, 0) == Approx(0.0).margin(1e-14));
    CIE_TEST_CHECK(vector(1, 0) == Approx(1.0));
    CIE_TEST_CHECK(vector(2, 0) == Approx(0.0).margin(1e-14));
}


CIE_TEST_CASE("Quaternion operator overloads", "[quaternions]")
{
    Quaternion<double> q({0.0,1.0,2.0,3.0});
    Quaternion<double> p({10.0,9.0,8.0,7.0});
    Quaternion<double> r = q + p;

    for (size_t i=0; i<4; ++i)
        CIE_TEST_CHECK(r[i] == Approx(10.0));

    r = q - p;
    for (size_t i=0; i<4; ++i)
        CIE_TEST_CHECK(r[i] == Approx(2*i-10.0));

    r *= 2.0;
    for (Size i=0; i<4; ++i)
        CIE_TEST_CHECK(r[i] == Approx(2.0 * (2*i-10.0)));

    r = p;
    r.conjugate();
    CIE_TEST_CHECK(r[0] == Approx(p[0]));
    for (size_t i=1; i<4; ++i)
        CIE_TEST_CHECK(r[i] == Approx(-p[i]));
}


} // namespace cie::linalg