// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/maths/inc/AffineTransform.hpp"


namespace cie::fem::maths {


CIE_TEST_CASE("AffineTransform", "[maths]")
{
    CIE_TEST_CASE_INIT("AffineTransform")

    {
        CIE_TEST_CASE_INIT("2D")
        //     (-1, 1)   3----|----2   ( 1, 1)
        //               |    |    |
        //          ----------o----------
        //               |    |    |
        //     (-1,-1)   0----|----1   ( 1,-1)
        //
        //               |||
        //               \ /
        //                v
        //
        // ( 2, 3)   2-----------------1   ( 4, 3)
        //   |       |                 |
        //   |       |                 |
        // ( 2, 1)   3-----------------0   ( 4, 1)
        //   |
        // --o--------------------------------
        //   |
        //   |
        using Transform = AffineTransform<double, 2>;
        using Point = Kernel<2,double>::Point;

        std::vector<Point> locals {
            {-1.0,-1.0},
            { 1.0,-1.0},
            {-1.0, 1.0},

            { 1.0, 1.0}
        };

        std::vector<Point> transformed {
            { 4.0, 1.0},
            { 4.0, 3.0},
            { 2.0, 1.0},

            { 2.0, 3.0}
        };

        Transform transform;
        {
            CIE_TEST_CASE_INIT("construction")
            CIE_TEST_CHECK_NOTHROW(transform = Transform(transformed.begin(), transformed.begin() + 3));
        }

        {
            CIE_TEST_CASE_INIT("transformation")

            for (Size i_point=0; i_point<locals.size(); ++i_point)
            {
                Point point;
                CIE_TEST_CHECK_NOTHROW(transform.evaluate(locals[i_point].begin(),
                                                          locals[i_point].end(),
                                                          point.begin()));
                CIE_TEST_REQUIRE(point.size() == 2);

                for (Size i_component=0; i_component<point.size(); ++i_component)
                    CIE_TEST_CHECK(point[i_component] == Approx(transformed[i_point][i_component]));
            }
        } // "transformation"

        {
            /*
            CIE_TEST_CASE_INIT("derivative")

            Transform::derivative_ptr p_derivative;
            CIE_TEST_REQUIRE_NOTHROW(p_derivative = transform.derivative());
            CIE_TEST_REQUIRE(p_derivative);

            Transform::derivative_type::value_type jacobian;
            jacobian(0,0) = 1.5;    jacobian(0,1) = 0.0;
            jacobian(1,0) = 0.0;    jacobian(1,1) = 1.5;

            for (const auto& r_point : transformed)
            {
                Transform::derivative_type::value_type test;
                CIE_TEST_CHECK_NOTHROW(test = p_derivative->evaluate(r_point));
                CIE_TEST_REQUIRE(test.rowSize() == 2);
                CIE_TEST_REQUIRE(test.columnSize() == 2);

                for (Size i_row=0; i_row<test.rowSize(); ++i_row)
                    for (Size i_column=0; i_column<test.columnSize(); i_column++)
                        CIE_TEST_CHECK(test(i_row, i_column) == Approx(jacobian(i_row, i_column)).margin(1e-14));
            }
            */
        } // "derivative"

        {
            CIE_TEST_CASE_INIT("inverse")

            const auto inverseTransform = transform.makeInverse();

            for (Size i_point=0; i_point<transformed.size(); ++i_point) {
                Point point;
                inverseTransform.evaluate(transformed[i_point].begin(),
                                          transformed[i_point].end(),
                                          point.begin());
                const auto& r_reference = locals[i_point];
                for (Size i_component=0; i_component<point.size(); ++i_component)
                    CIE_TEST_CHECK(point[i_component] == Approx(r_reference[i_component]));
            }
        } // "inverse"
    } // "2D"


    {
        CIE_TEST_CASE_INIT("3D")
        using Transform = AffineTransform<double, 3>;
        using Point = Kernel<3,double>::Point;

        std::vector<Point> locals {
            {-1, -1, -1},
            {1, -1, -1},
            {-1, 1, -1},
            {-1, -1, 1},

            {1, 1, 1}
        };

        std::vector<Point> transformed {
            {-1, -1, -1},
            {2, -1, -1},
            {-1, 4, -1},
            {-1, -1, 6},

            {2, 4, 6}
        };

        Transform transform;
        {
            CIE_TEST_CASE_INIT("construction")
            CIE_TEST_CHECK_NOTHROW(transform = Transform(transformed.begin(), transformed.begin() + 4));
        }

        {
            CIE_TEST_CASE_INIT("transformation")

            for (Size i_point=0; i_point<locals.size(); ++i_point) {
                Point point;
                CIE_TEST_CHECK_NOTHROW(transform.evaluate(locals[i_point].begin(),
                                                          locals[i_point].end(),
                                                          point.begin()));
                CIE_TEST_REQUIRE(point.size() == 3);
                for (Size i_component=0; i_component<point.size(); ++i_component)
                    CIE_TEST_CHECK(point[i_component] == Approx(transformed[i_point][i_component]));
            }
        } // "transformation"

        {
            /*
            CIE_TEST_CASE_INIT("derivative")

            Transform::derivative_ptr p_derivative;
            CIE_TEST_REQUIRE_NOTHROW(p_derivative = transform.derivative());
            CIE_TEST_REQUIRE(p_derivative);

            Transform::derivative_type::value_type jacobian;
            jacobian(0,0) = 1.5;    jacobian(0,1) = 0.0;    jacobian(0,2) = 0.0;
            jacobian(1,0) = 0.0;    jacobian(1,1) = 2.5;    jacobian(1,2) = 0.0;
            jacobian(2,0) = 0.0;    jacobian(2,1) = 0.0;    jacobian(2,2) = 3.5;

            for (const auto& r_point : transformed)
            {
                Transform::derivative_type::value_type test;
                CIE_TEST_CHECK_NOTHROW(test = p_derivative->evaluate(r_point));
                CIE_TEST_REQUIRE(test.rowSize() == 3);
                CIE_TEST_REQUIRE(test.columnSize() == 3);

                for (Size i_row=0; i_row<test.rowSize(); ++i_row)
                    for (Size i_column=0; i_column<test.columnSize(); i_column++)
                        CIE_TEST_CHECK(test(i_row, i_column) == Approx(jacobian(i_row, i_column)).margin(1e-14));
            }
            */
        } // "derivative"

        {
            CIE_TEST_CASE_INIT("inverse")

            const auto inverseTransform = transform.makeInverse();

            for (Size i_point=0; i_point<transformed.size(); ++i_point) {
                Point point;
                inverseTransform.evaluate(transformed[i_point].begin(),
                                          transformed[i_point].end(),
                                          point.begin());
                const auto& r_reference = locals[i_point];
                for (Size i_component=0; i_component<point.size(); ++i_component)
                    CIE_TEST_CHECK(point[i_component] == Approx(r_reference[i_component]));
            }
        } // "inverse"
    } // "3D"
}


} // namespace cie::fem::maths
