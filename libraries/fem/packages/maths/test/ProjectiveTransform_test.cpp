// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/ProjectiveTransform.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- STL Includes ---
//#include <iostream>


namespace cie::fem::maths {


CIE_TEST_CASE("ProjectiveTransform", "[maths]")
{
    CIE_TEST_CASE_INIT("ProjectiveTransform")

    {
        CIE_TEST_CASE_INIT("2D")

        constexpr unsigned Dimension = 2u;
        using Point = Kernel<Dimension,double>::Point;
        using Transform = ProjectiveTransform<double,Dimension>;
        CIE_TEST_CHECK(SpatialTransform<Transform>);

        StaticArray<Point,4> transformedPoints {{1.0, 1.0},
                                                {3.0, 3.0},
                                                {0.0, 1.0},
                                                {3.0, 4.0}};
        const auto transform = Transform(transformedPoints.begin(),
                                         transformedPoints.end());
        const auto jacobian = transform.makeDerivative();


        const double delta = 1e-8;

        for (const auto& input : {Point {-1.0, -1.0},
                                  Point { 1.0, -1.0},
                                  Point {-1.0,  1.0},
                                  Point { 1.0,  1.0},
                                  Point { 0.0, -1.0},
                                  Point { 1.0, -1.0},
                                  Point { 0.0,  0.0},
                                  Point { 1.0,  0.0},
                                  Point {-1.2,  3.4},
                                  Point { 3.4,  3.4},
                                  Point {-1.2, -1.2},
                                  Point { 3.4, -1.2}}) {
            Point inputDelta;
            Eigen::Matrix<double,Dimension,Dimension,Eigen::RowMajor> output;
            Eigen::Matrix<double,Dimension,Dimension> outputBase, outputDelta;

            jacobian.evaluate(input.data(),
                            input.data() + 2,
                            output.data());

            inputDelta = input;
            inputDelta[0] += delta;
            transform.evaluate(input.data(),
                               input.data() + 2,
                               outputBase.data());
            transform.evaluate(inputDelta.begin(),
                               inputDelta.end(),
                               outputDelta.data());

            inputDelta = input;
            inputDelta[1] += delta;
            transform.evaluate(input.data(),
                               input.data() + 2,
                               outputBase.data() + 2);
            transform.evaluate(inputDelta.begin(),
                               inputDelta.end(),
                               outputDelta.data() + 2);

            // @todo incorrect derivative (impl does what it should, the derivation is wrong)
            const Eigen::Matrix<double,2,2> reference = ((outputDelta - outputBase) / delta);
            for (unsigned i_row=0; i_row<2; ++i_row) {
                for (unsigned i_column=0; i_column<2; ++i_column) {
                    CIE_TEST_CHECK(output(i_row, i_column) == Approx(reference(i_row, i_column)).margin(1e-5));
                }
            }
        }
    }

    {
        CIE_TEST_CASE_INIT("3D")

        constexpr unsigned Dimension = 3u;
        using Point = Kernel<Dimension,double>::Point;
        using Transform = ProjectiveTransform<double,Dimension>;
        CIE_TEST_CHECK(SpatialTransform<Transform>);

        StaticArray<Point,8> transformedPoints {{-1.0, -1.0, -1.0},
                                                { 1.0, -1.0, -1.0},
                                                {-1.0,  1.0, -1.0},
                                                { 1.0,  1.0, -1.0},
                                                {-1.0, -1.0,  1.0},
                                                { 1.0, -1.0,  1.0},
                                                {-1.0,  1.0,  1.0},
                                                { 1.0,  1.0,  1.0}};
        const auto transform = Transform(transformedPoints.begin(),
                                         transformedPoints.end());
        [[maybe_unused]] const auto jacobian = transform.makeDerivative();


        const double delta = 1e-8;

        for (const auto& input : {Point {-1.0, -1.0, -1.0},
                                  Point { 1.0, -1.0, -1.0},
                                  Point {-1.0,  1.0, -1.0},
                                  Point { 1.0,  1.0, -1.0},
                                  Point {-1.0, -1.0,  1.0},
                                  Point { 1.0, -1.0,  1.0},
                                  Point {-1.0,  1.0,  1.0},
                                  Point { 1.0,  1.0,  1.0}}) {
            Point inputDelta;
            Eigen::Matrix<double,Dimension,Dimension,Eigen::RowMajor> output;
            Eigen::Matrix<double,Dimension,Dimension> outputBase, outputDelta;

            jacobian.evaluate(input.data(),
                              input.data() + input.size(),
                              output.data());

            for (unsigned iDimension=0; iDimension<Dimension; ++iDimension) {
                CIE_TEST_CHECK_NOTHROW(transform.evaluate(
                    input.data(),
                    input.data(),
                    outputBase.data() + iDimension * Dimension
                ));

                inputDelta = input;
                inputDelta[iDimension] += delta;

                CIE_TEST_CHECK_NOTHROW(transform.evaluate(
                    inputDelta.data(),
                    inputDelta.data(),
                    outputDelta.data() + iDimension * Dimension
                ));
            } // for iDimension in range(Dimension)

            // @todo incorrect derivative (impl does what it should, the derivation is wrong)
            const Eigen::Matrix<double,Dimension,Dimension> reference = ((outputDelta - outputBase) / delta);
            for (unsigned i_row=0; i_row<Dimension; ++i_row) {
                for (unsigned i_column=0; i_column<Dimension; ++i_column) {
                    CIE_TEST_CHECK(output(i_row, i_column) == Approx(reference(i_row, i_column)).margin(1e-5));
                }
            }
        }
    }
}


} // namespace cie::fem::maths
