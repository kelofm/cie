// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/AffineEmbedding.hpp"

// --- STL Includes ---
#include <numeric>


namespace cie::fem::maths {


CIE_TEST_CASE("AffineEmbedding", "[maths]")
{
    CIE_TEST_CASE_INIT("AffineEmbedding")

    {
        CIE_TEST_CASE_INIT("1D => 2D")
        using Embedding = AffineEmbedding<double,1,2>;

        StaticArray<double,1> in;
        StaticArray<double,2> out;

        // Check default-constructed embedding.
        CIE_TEST_CHECK_NOTHROW(Embedding());
        Embedding embedding;

        in[0] = -1.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(-1.0));
        CIE_TEST_CHECK(out[1] == Approx(-1.0));

        in[0] = 0.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(0.0));
        CIE_TEST_CHECK(out[1] == Approx(-1.0));

        in[0] = 1.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(1.0));
        CIE_TEST_CHECK(out[1] == Approx(-1.0));

        // Check embedding from 2 endpoints.
        StaticArray<Embedding::OutPoint,2> transformed;
        transformed[0] = Embedding::OutPoint {10.0, 5.0};
        transformed[1] = Embedding::OutPoint {-2.0, 15.0};
        CIE_TEST_CHECK_NOTHROW(embedding = Embedding(transformed));

        in[0] = -1.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(10.0));
        CIE_TEST_CHECK(out[1] == Approx(5.0));

        in[0] = 0.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(4.0));
        CIE_TEST_CHECK(out[1] == Approx(10.0));

        in[0] = 1.0;
        CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
        CIE_TEST_CHECK(out[0] == Approx(-2.0));
        CIE_TEST_CHECK(out[1] == Approx(15.0));

        // Check the embedding's derivative.
        {
            CIE_TEST_REQUIRE_NOTHROW(embedding.makeDerivative());
            const auto jacobian = embedding.makeDerivative();

            const StaticArray<double,2> segment {
                transformed[1][0] - transformed[0][0],
                transformed[1][1] - transformed[0][1]
            };
            const double segmentNorm = std::sqrt(std::inner_product(
                segment.begin(),
                segment.end(),
                segment.begin(),
                0.0));

            CIE_TEST_CHECK_NOTHROW(embedding.evaluate(in, out));
            StaticArray<double,1> reference, perturbed;
            constexpr double perturbationNorm = 1e0;
            StaticArray<double,2> transformedPerturbed;

            perturbed = in;
            perturbed[0] += perturbationNorm;
            CIE_TEST_CHECK_NOTHROW(embedding.evaluate(perturbed, transformedPerturbed));
            reference[0] = (transformedPerturbed[0] - out[0]) / perturbationNorm;
            reference[1] = (transformedPerturbed[1] - out[1]) / perturbationNorm;

            StaticArray<double,2> test;
            CIE_TEST_CHECK_NOTHROW(jacobian.evaluate(in, test));
            CIE_TEST_CHECK(test[0] == Approx(reference[0]));
            CIE_TEST_CHECK(test[1] == Approx(reference[1]));

            CIE_TEST_CHECK(jacobian.evaluateDeterminant(out) == Approx(segmentNorm / 2.0));
        }

        // Check the embedding's inverse.
        CIE_TEST_REQUIRE_NOTHROW(embedding.makeInverse());
        const auto projection = embedding.makeInverse();

        out = StaticArray<double,2> { 10.0,  5.0};
        CIE_TEST_CHECK_NOTHROW(projection.evaluate(out, in));
        CIE_TEST_CHECK(in[0] == Approx(-1.0));

        out = StaticArray<double,2> {  4.0, 10.0};
        CIE_TEST_CHECK_NOTHROW(projection.evaluate(out, in));
        CIE_TEST_CHECK(in[0] == Approx(0.0).margin(1e-14));

        out = StaticArray<double,2> { -2.0, 15.0};
        CIE_TEST_CHECK_NOTHROW(projection.evaluate(out, in));
        CIE_TEST_CHECK(in[0] == Approx(1.0));

        // Check the embedding's inverse when the input point
        // does not lie on the transformed line segment.
        out = StaticArray<double,2> {  0.0,  0.0};
        CIE_TEST_CHECK_NOTHROW(projection.evaluate(out, in));
        CIE_TEST_CHECK(in[0] == Approx(-0.42622950819672134));

        // Check the derivative of the embedding's inverse.
        {
            CIE_TEST_CHECK_NOTHROW(projection.makeDerivative());
            const auto jacobian = projection.makeDerivative();

            const StaticArray<double,2> segment {
                transformed[1][0] - transformed[0][0],
                transformed[1][1] - transformed[0][1]
            };
            const double segmentNorm = std::sqrt(std::inner_product(segment.begin(),
                                                                    segment.end(),
                                                                    segment.begin(),
                                                                    0.0));

            CIE_TEST_CHECK_NOTHROW(projection.evaluate(out, in));
            StaticArray<double,2> reference, perturbed;
            constexpr double perturbationNorm = 1e0;
            StaticArray<double,1> transformedPerturbed;

            perturbed = out;
            perturbed[0] += perturbationNorm;
            CIE_TEST_CHECK_NOTHROW(projection.evaluate(perturbed, transformedPerturbed));
            reference[0] = (transformedPerturbed[0] - in[0]) / perturbationNorm;

            perturbed = out;
            perturbed[1] += perturbationNorm;
            CIE_TEST_CHECK_NOTHROW(projection.evaluate(perturbed, transformedPerturbed));
            reference[1] = (transformedPerturbed[0] - in[0]) / perturbationNorm;

            StaticArray<double,2> test;
            CIE_TEST_CHECK_NOTHROW(jacobian.evaluate(out, test));
            CIE_TEST_CHECK(test[0] == Approx(reference[0]));
            CIE_TEST_CHECK(test[1] == Approx(reference[1]));

            CIE_TEST_CHECK(jacobian.evaluateDeterminant(out) == Approx(2.0 / segmentNorm));
        }
    }
}


} // namespace cie::fem::maths
