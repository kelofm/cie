// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/kernel.hpp"

namespace cie::fem::maths {


CIE_TEST_CASE( "AnsatzSpace", "[maths]" ) {
    CIE_TEST_CASE_INIT( "AnsatzSpace" )

    {
        CIE_TEST_CASE_INIT("dynamic")
        using Point = Kernel<2,double>::Point;

        #define CIE_TMP_CHECK(ANSATZ_SPACE)                                             \
            {                                                                           \
                double x = 11.0;                                                        \
                double y = 9.0;                                                         \
                Point p0 {-x, -y};                                                      \
                Point p1 {-x,  y};                                                      \
                Point p2 { x, -y};                                                      \
                Point p3 { x,  y};                                                      \
                StaticArray<double,9> results;                                          \
                                                                                        \
                CIE_TEST_REQUIRE(ANSATZ_SPACE.size() == 9);                             \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p0, results));             \
                CIE_TEST_CHECK(results[0] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[1] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[2] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[3] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[4] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[5] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[6] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[7] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p1, results));             \
                CIE_TEST_CHECK(results[0] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[1] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[2] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[3] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[4] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[5] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[6] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[7] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p2, results));             \
                CIE_TEST_CHECK(results[0] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[1] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[2] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[3] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[4] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[5] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[6] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[7] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p3, results));             \
                CIE_TEST_CHECK(results[0] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[1] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[2] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[3] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[4] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[5] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[6] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[7] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
            }

        using Basis = Polynomial<double>;
        using Ansatz = AnsatzSpace<Basis,2>;

        // Check default constructor.
        CIE_TEST_CHECK_NOTHROW(Ansatz());

        const DynamicArray<Basis> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5}),
            Basis(Basis::Coefficients { 0.5,  0.5}),
            Basis(Basis::Coefficients { 0.0,  0.0, 1.0})};

        // Check constructor.
        {
            CIE_TEST_REQUIRE_NOTHROW(Ansatz(basisFunctions));
            Ansatz ansatzSpace(basisFunctions);
            CIE_TMP_CHECK(ansatzSpace)
        }

        // Check move constructor.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check copy constructor.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other = *pSwap;
            CIE_TMP_CHECK((*pSwap))
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check move assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other;
            other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check copy assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other;
            other = *pSwap;
            CIE_TMP_CHECK((*pSwap))
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        #undef CIE_TMP_CHECK
    }

    {
        CIE_TEST_CASE_INIT("static")
        using Point = Kernel<2,double>::Point;

        #define CIE_TMP_CHECK(ANSATZ_SPACE)                                             \
            {                                                                           \
                double x = 11.0;                                                        \
                double y = 9.0;                                                         \
                Point p0 {-x, -y};                                                      \
                Point p1 {-x,  y};                                                      \
                Point p2 { x, -y};                                                      \
                Point p3 { x,  y};                                                      \
                StaticArray<double,9> results;                                          \
                                                                                        \
                CIE_TEST_REQUIRE(ANSATZ_SPACE.size() == 9);                             \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p0, results));             \
                CIE_TEST_CHECK(results[0] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[1] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[2] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[3] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[4] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[5] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[6] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[7] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p1, results));             \
                CIE_TEST_CHECK(results[0] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[1] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[2] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[3] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[4] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[5] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[6] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[7] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p2, results));             \
                CIE_TEST_CHECK(results[0] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[1] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[2] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[3] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[4] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[5] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[6] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[7] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p3, results));             \
                CIE_TEST_CHECK(results[0] == Approx(  20.0));                           \
                CIE_TEST_CHECK(results[1] == Approx( -24.0));                           \
                CIE_TEST_CHECK(results[2] == Approx(-484.0));                           \
                CIE_TEST_CHECK(results[3] == Approx( -25.0));                           \
                CIE_TEST_CHECK(results[4] == Approx(  30.0));                           \
                CIE_TEST_CHECK(results[5] == Approx( 605.0));                           \
                CIE_TEST_CHECK(results[6] == Approx(-405.0));                           \
                CIE_TEST_CHECK(results[7] == Approx( 486.0));                           \
                CIE_TEST_CHECK(results[8] == Approx(9801.0));                           \
            }

        using Basis = Polynomial<double,2>;
        using Ansatz = AnsatzSpace<Basis,2,3>;

        // Check default constructor.
        CIE_TEST_CHECK_NOTHROW(Ansatz());

        const std::array<Basis,3> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5,  0.0}),
            Basis(Basis::Coefficients { 0.5,  0.5,  0.0}),
            Basis(Basis::Coefficients { 0.0,  0.0,  1.0})};

        // Check constructor.
        {
            CIE_TEST_REQUIRE_NOTHROW(Ansatz(basisFunctions));
            Ansatz ansatzSpace(basisFunctions);
            CIE_TMP_CHECK(ansatzSpace)
        }

        // Check move constructor.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check copy constructor.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other = *pSwap;
            CIE_TMP_CHECK((*pSwap))
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check move assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other;
            other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        // Check copy assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz>(basisFunctions);
            Ansatz other;
            other = *pSwap;
            CIE_TMP_CHECK((*pSwap))
            pSwap.reset();
            CIE_TMP_CHECK(other)
        }

        #undef CIE_TMP_CHECK
    }
}


CIE_TEST_CASE( "AnsatzSpaceDerivative", "[maths]" ) {
    CIE_TEST_CASE_INIT( "AnsatzSpaceDerivative" )

    {
        CIE_TEST_CASE_INIT("dynamic")
        using Basis = Polynomial<double>;
        using Ansatz = AnsatzSpace<Basis,2>;

        const DynamicArray<Basis> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5}),
            Basis(Basis::Coefficients { 0.5,  0.5})};

        #define CIE_TMP_CHECK(ANSATZ_DERIVATIVE)                                                \
            {                                                                                   \
                CIE_TEST_REQUIRE(ANSATZ_DERIVATIVE.size() == 8);                                \
                                                                                                \
                using Point = Kernel<2,double>::Point;                                          \
                const double x = 2;                                                             \
                const double y = 3;                                                             \
                const Point p0 {-x, -y};                                                        \
                const Point p1 { x, -y};                                                        \
                const Point p2 {-x,  y};                                                        \
                const Point p3 { x,  y};                                                        \
                                                                                                \
                Eigen::Matrix<double,4,2> buffer;                                               \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p0, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p1, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p2, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p3, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                  \
            }

        // Check default constructor.
        CIE_TEST_REQUIRE_NOTHROW(Ansatz::Derivative());

        // Check constructor.
        {
            CIE_TEST_REQUIRE_NOTHROW(Ansatz::Derivative(basisFunctions));
            const auto derivative = Ansatz::Derivative(basisFunctions);
            CIE_TMP_CHECK(derivative)
        }

        // Check move constructor.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check copy constructor.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other = *pSwap;
            CIE_TMP_CHECK((*pSwap));
            CIE_TMP_CHECK(other);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check move assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other;
            other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check copy assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other;
            other = *pSwap;
            CIE_TMP_CHECK((*pSwap));
            CIE_TMP_CHECK(other);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }
    }

    {
        CIE_TEST_CASE_INIT("static")
        using Basis = Polynomial<double,1>;
        using Ansatz = AnsatzSpace<Basis,2,2>;

        const std::array<Basis,2> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5}),
            Basis(Basis::Coefficients { 0.5,  0.5})};

        #define CIE_TMP_CHECK(ANSATZ_DERIVATIVE) \
            {                                                                                   \
                CIE_TEST_REQUIRE(ANSATZ_DERIVATIVE.size() == 8);                                \
                                                                                                \
                using Point = Kernel<2,double>::Point;                                          \
                const double x = 2;                                                             \
                const double y = 3;                                                             \
                const Point p0 {-x, -y};                                                        \
                const Point p1 { x, -y};                                                        \
                const Point p2 {-x,  y};                                                        \
                const Point p3 { x,  y};                                                        \
                                                                                                \
                Eigen::Matrix<double,4,2> buffer;                                               \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p0, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p1, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p2, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                  \
                                                                                                \
                ANSATZ_DERIVATIVE.evaluate(p3, {buffer.data(), buffer.data() + buffer.size()}); \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                  \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                  \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                  \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                  \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                  \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                  \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                  \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                  \
            }

        // Check default constructor.
        CIE_TEST_REQUIRE_NOTHROW(Ansatz::Derivative());

        // Check constructor.
        {
            CIE_TEST_REQUIRE_NOTHROW(Ansatz::Derivative(basisFunctions));
            const Ansatz::Derivative derivative(basisFunctions);
            CIE_TMP_CHECK(derivative)
        }

        // Check move constructor.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check copy constructor.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other = *pSwap;
            CIE_TMP_CHECK((*pSwap));
            CIE_TMP_CHECK(other);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check move assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other;
            other = std::move(*pSwap);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }

        // Check copy assignment operator.
        {
            auto pSwap = std::make_unique<Ansatz::Derivative>(Ansatz::Derivative(basisFunctions));
            Ansatz::Derivative other;
            other = *pSwap;
            CIE_TMP_CHECK((*pSwap));
            CIE_TMP_CHECK(other);
            pSwap.reset();
            CIE_TMP_CHECK(other);
        }
    }
}


} // namespace cie::fem::maths
