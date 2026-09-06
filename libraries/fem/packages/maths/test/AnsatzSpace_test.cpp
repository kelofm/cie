// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- STL Includes ---
#include <sstream>


namespace cie::fem::maths {


CIE_TEST_CASE( "AnsatzSpace", "[maths]" ) {
    CIE_TEST_CASE_INIT( "AnsatzSpace" )

    using T = double;
    std::allocator<T> allocator;

    {
        CIE_TEST_CASE_INIT("dynamic")
        using Point = Kernel<2,T>::Point;
        using Basis = Polynomial<double>;
        using Ansatz = AnsatzSpace<Basis,2>;

        #define CIE_TMP_CHECK(ANSATZ_SPACE)                                             \
            {                                                                           \
                double x = 11.0;                                                        \
                double y = 9.0;                                                         \
                Point p0 {-x, -y};                                                      \
                Point p1 {-x,  y};                                                      \
                Point p2 { x, -y};                                                      \
                Point p3 { x,  y};                                                      \
                StaticArray<double,9> results;                                          \
                std::vector<double> buffer(ANSATZ_SPACE.bufferSize());                  \
                                                                                        \
                CIE_TEST_REQUIRE(ANSATZ_SPACE.size() == 9);                             \
                                                                                        \
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p0, results, buffer));     \
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
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p1, results, buffer));     \
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
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p2, results, buffer));     \
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
                CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(p3, results, buffer));     \
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

        // Check default constructor.
        CIE_TEST_CHECK_NOTHROW(Ansatz());

        const Ansatz::AnsatzSet basisFunctions {
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

        // Check serialization and deserialization.
        {
            Ansatz instance(basisFunctions);
            std::stringstream stream;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::serialize(
                stream,
                instance));
            Ansatz other;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::deserialize(
                stream,
                other,
                allocator));
            CIE_TMP_CHECK(other);
            CIE_TMP_CHECK(instance);
        }

        #undef CIE_TMP_CHECK
    }

    {
        CIE_TEST_CASE_INIT("static")
        using Point = Kernel<2,T>::Point;

        #define CIE_TMP_CHECK(ANSATZ_SPACE)                                                                         \
            {                                                                                                       \
                T x = 11.0;                                                                                         \
                T y = 9.0;                                                                                          \
                StaticArray<T,9> results;                                                                           \
                std::vector<T> buffer(ANSATZ_SPACE.bufferSize());                                                   \
                                                                                                                    \
                Ansatz::Context slice;                                                                              \
                slice.data()[0] = 0;                                                                                \
                slice.data()[1] = 2;                                                                                \
                slice.data()[2] = 2;                                                                                \
                slice.data()[3] = 3;                                                                                \
                std::array<Ansatz::Context::Value,2> sliceIndices;                                                  \
                std::array<T,2> slicedResults;                                                                      \
                const std::array<Ansatz::Context::Value,2> range {3, 3};                                            \
                CIE_TEST_REQUIRE_NOTHROW(slice.expand(range, sliceIndices));                                        \
                CIE_TEST_REQUIRE(sliceIndices[0] == 6);                                                             \
                CIE_TEST_REQUIRE(sliceIndices[1] == 7);                                                             \
                CIE_TEST_REQUIRE(ANSATZ_SPACE.size() == 9);                                                         \
                CIE_TEST_REQUIRE(ANSATZ_SPACE.size(slice) == 2);                                                    \
                                                                                                                    \
                {                                                                                                   \
                    Point point {-x, -y};                                                                           \
                    const std::array<T,9> reference {                                                               \
                          30.0,  -25.0,  605.0,                                                                     \
                         -24.0,   20.0, -484.0,                                                                     \
                         486.0, -405.0, 9801.0};                                                                    \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, results, buffer));                                                                   \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, slicedResults, buffer, slice));                                                      \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                                                                                                                    \
                {                                                                                                   \
                    Point point {-x,  y};                                                                           \
                    const std::array<T,9> reference {                                                               \
                         -24.0,   20.0, -484.0,                                                                     \
                          30.0,  -25.0,  605.0,                                                                     \
                         486.0, -405.0, 9801.0};                                                                    \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, results, buffer));                                                                   \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, slicedResults, buffer, slice));                                                      \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                                                                                                                    \
                {                                                                                                   \
                    Point point { x, -y};                                                                           \
                    const std::array<T,9> reference {                                                               \
                         -25.0,   30.0,  605.0,                                                                     \
                          20.0,  -24.0, -484.0,                                                                     \
                        -405.0,  486.0, 9801.0};                                                                    \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, results, buffer));                                                                   \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, slicedResults, buffer, slice));                                                      \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                                                                                                                    \
                {                                                                                                   \
                    Point point { x,  y};                                                                           \
                    const std::array<T,9> reference {                                                               \
                          20.0,  -24.0, -484.0,                                                                     \
                         -25.0,   30.0,  605.0,                                                                     \
                        -405.0,  486.0, 9801.0};                                                                    \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, results, buffer));                                                                   \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_SPACE.evaluate(                                                   \
                        point, slicedResults, buffer, slice));                                                      \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
            }

        using Basis = Polynomial<T,2>;
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

        // Check serialization and deserialization.
        {
            Ansatz instance(basisFunctions);
            std::stringstream stream;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::serialize(
                stream,
                instance));
            Ansatz other;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::deserialize(
                stream,
                other,
                allocator));
            CIE_TMP_CHECK(other);
            CIE_TMP_CHECK(instance);
        }

        #undef CIE_TMP_CHECK
    }
}


CIE_TEST_CASE( "AnsatzSpaceDerivative", "[maths]" ) {
    CIE_TEST_CASE_INIT( "AnsatzSpaceDerivative" )

    using T = double;
    std::allocator<T> allocator;

    {
        CIE_TEST_CASE_INIT("dynamic")
        using Basis = Polynomial<T>;
        using Ansatz = AnsatzSpace<Basis,2>;

        const DynamicArray<Basis> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5}),
            Basis(Basis::Coefficients { 0.5,  0.5})};

        #define CIE_TMP_CHECK(ANSATZ_DERIVATIVE)                                                    \
            {                                                                                       \
                CIE_TEST_REQUIRE(ANSATZ_DERIVATIVE.size() == 8);                                    \
                                                                                                    \
                using Point = Kernel<2,T>::Point;                                                   \
                const T x = 2;                                                                      \
                const T y = 3;                                                                      \
                const Point p0 {-x, -y};                                                            \
                const Point p1 { x, -y};                                                            \
                const Point p2 {-x,  y};                                                            \
                const Point p3 { x,  y};                                                            \
                std::vector<T> b(ANSATZ_DERIVATIVE.bufferSize());                                   \
                                                                                                    \
                Eigen::Matrix<T,4,2> buffer;                                                        \
                                                                                                    \
                ANSATZ_DERIVATIVE.evaluate(p0, {buffer.data(), buffer.data() + buffer.size()}, b);  \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                      \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                      \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                      \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                      \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                      \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                      \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                      \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                      \
                                                                                                    \
                ANSATZ_DERIVATIVE.evaluate(p1, {buffer.data(), buffer.data() + buffer.size()}, b);  \
                CIE_TEST_CHECK(buffer(0, 0) == Approx(-1.00));                                      \
                CIE_TEST_CHECK(buffer(1, 0) == Approx( 1.00));                                      \
                CIE_TEST_CHECK(buffer(2, 0) == Approx( 0.50));                                      \
                CIE_TEST_CHECK(buffer(3, 0) == Approx(-0.50));                                      \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                      \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                      \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                      \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                      \
                                                                                                    \
                ANSATZ_DERIVATIVE.evaluate(p2, {buffer.data(), buffer.data() + buffer.size()}, b);  \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                      \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                      \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                      \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                      \
                CIE_TEST_CHECK(buffer(0, 1) == Approx(-0.75));                                      \
                CIE_TEST_CHECK(buffer(1, 1) == Approx( 0.25));                                      \
                CIE_TEST_CHECK(buffer(2, 1) == Approx( 0.75));                                      \
                CIE_TEST_CHECK(buffer(3, 1) == Approx(-0.25));                                      \
                                                                                                    \
                ANSATZ_DERIVATIVE.evaluate(p3, {buffer.data(), buffer.data() + buffer.size()}, b);  \
                CIE_TEST_CHECK(buffer(0, 0) == Approx( 0.50));                                      \
                CIE_TEST_CHECK(buffer(1, 0) == Approx(-0.50));                                      \
                CIE_TEST_CHECK(buffer(2, 0) == Approx(-1.00));                                      \
                CIE_TEST_CHECK(buffer(3, 0) == Approx( 1.00));                                      \
                CIE_TEST_CHECK(buffer(0, 1) == Approx( 0.25));                                      \
                CIE_TEST_CHECK(buffer(1, 1) == Approx(-0.75));                                      \
                CIE_TEST_CHECK(buffer(2, 1) == Approx(-0.25));                                      \
                CIE_TEST_CHECK(buffer(3, 1) == Approx( 0.75));                                      \
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

        // Check serialization and deserialization.
        {
            Ansatz::Derivative instance(basisFunctions);
            std::stringstream stream;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::serialize(
                stream,
                instance));
            Ansatz::Derivative other;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::deserialize(
                stream,
                other,
                allocator));
            CIE_TMP_CHECK(other);
            CIE_TMP_CHECK(instance);
        }

        #undef CIE_TMP_CHECK
    }

    {
        CIE_TEST_CASE_INIT("static")
        using Basis = Polynomial<T,1>;
        using Ansatz = AnsatzSpace<Basis,2,2>;

        const std::array<Basis,2> basisFunctions {
            Basis(Basis::Coefficients { 0.5, -0.5}),
            Basis(Basis::Coefficients { 0.5,  0.5})};

        #define CIE_TMP_CHECK(ANSATZ_DERIVATIVE) \
            {                                                                                                       \
                CIE_TEST_REQUIRE(ANSATZ_DERIVATIVE.size() == 8);                                                    \
                                                                                                                    \
                using Point = Kernel<2,T>::Point;                                                                   \
                const T x = 2;                                                                                      \
                const T y = 3;                                                                                      \
                const Point p0 {-x, -y};                                                                            \
                const Point p1 { x, -y};                                                                            \
                const Point p2 {-x,  y};                                                                            \
                const Point p3 { x,  y};                                                                            \
                std::vector<T> b(ANSATZ_DERIVATIVE.bufferSize());                                                   \
                                                                                                                    \
                std::array<T,8> results;                                                                            \
                Ansatz::Context slice;                                                                              \
                slice.data()[0] = 0;                                                                                \
                slice.data()[1] = 1;                                                                                \
                slice.data()[2] = 1;                                                                                \
                slice.data()[3] = 2;                                                                                \
                std::array<Ansatz::Context::Value,2> sliceIndices;                                                  \
                std::array<T,2> slicedResults;                                                                      \
                const std::array<Ansatz::Context::Value,2> range {2, 2};                                            \
                CIE_TEST_REQUIRE_NOTHROW(slice.expand(range, sliceIndices));                                        \
                CIE_TEST_REQUIRE(sliceIndices[0] == 2);                                                             \
                sliceIndices[1] = (sliceIndices[0] + 1) * 2;                                                        \
                CIE_TEST_REQUIRE(ANSATZ_DERIVATIVE.size(slice) == 2);                                               \
                                                                                                                    \
                {                                                                                                   \
                    Point point {-x, -y};                                                                           \
                    const std::array<T,8> reference {                                                               \
                        -1.00,  1.00,  0.50, -0.50,                                                                 \
                        -0.75,  0.25,  0.75, -0.25};                                                                \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, results, b));                                                                        \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, slicedResults, b, slice));                                                           \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                {                                                                                                   \
                    Point point { x, -y};                                                                           \
                    const std::array<T,8> reference {                                                               \
                        -1.00,  1.00,  0.50, -0.50,                                                                 \
                         0.25, -0.75, -0.25,  0.75};                                                                \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, results, b));                                                                        \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, slicedResults, b, slice));                                                           \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                {                                                                                                   \
                    Point point {-x,  y};                                                                           \
                    const std::array<T,8> reference {                                                               \
                         0.50, -0.50, -1.00,  1.00,                                                                 \
                        -0.75,  0.25,  0.75, -0.25};                                                                \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, results, b));                                                                        \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, slicedResults, b, slice));                                                           \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
                {                                                                                                   \
                    Point point { x,  y};                                                                           \
                    const std::array<T,8> reference {                                                               \
                         0.50, -0.50, -1.00,  1.00,                                                                 \
                         0.25, -0.75, -0.25,  0.75};                                                                \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, results, b));                                                                        \
                    for (std::size_t iComponent=0ul; iComponent<reference.size(); ++iComponent)                     \
                        CIE_TEST_CHECK(results[iComponent] == Approx(reference[iComponent]));                       \
                                                                                                                    \
                    CIE_TEST_CHECK_NOTHROW(ANSATZ_DERIVATIVE.evaluate(                                              \
                        point, slicedResults, b, slice));                                                           \
                    for (std::size_t iComponent=0ul; iComponent<sliceIndices.size(); ++iComponent)                  \
                        CIE_TEST_CHECK(slicedResults[iComponent] == Approx(reference[sliceIndices[iComponent]]));   \
                }                                                                                                   \
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

        // Check serialization and deserialization.
        {
            Ansatz::Derivative instance(basisFunctions);
            std::stringstream stream;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::serialize(
                stream,
                instance));
            Ansatz::Derivative other;
            CIE_TEST_REQUIRE_NOTHROW(cie::io::Serializer<tags::Binary>::deserialize(
                stream,
                other,
                allocator));
            CIE_TMP_CHECK(other);
            CIE_TMP_CHECK(instance);
        }

        #undef CIE_TMP_CHECK
    }
}


} // namespace cie::fem::maths
