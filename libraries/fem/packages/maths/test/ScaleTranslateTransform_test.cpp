// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/utilities/inc/kernel.hpp"


namespace cie::fem::maths {


#define CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST(TransformType)                                                        \
    CIE_TEST_CASE(#TransformType, "[maths]")                                                                                \
    {                                                                                                                       \
        CIE_TEST_CASE_INIT(#TransformType)                                                                                  \
        /*                                                                                                                  \
        //     (-1, 1)   3----|----2   ( 1, 1)                                                                              \
        //               |    |    |                                                                                        \
        //          ----------o----------                                                                                   \
        //               |    |    |                                                                                        \
        //     (-1,-1)   0----|----1   ( 1,-1)                                                                              \
        //                                                                                                                  \
        //                   |||                                                                                            \
        //                   \ /                                                                                            \
        //                    v                                                                                             \
        //                                                                                                                  \
        //   | ( 2, 3)   2---------------3   ( 6, 3)                                                                        \
        //   |           |               |                                                                                  \
        //   |           |               |                                                                                  \
        //   |           |               |                                                                                  \
        //   | ( 2, 1)   1---------------0   ( 6, 1)                                                                        \
        //   |                                                                                                              \
        // --o------------------------                                                                                      \
        //   |                                                                                                              \
        //   |                                                                                                              \
        */                                                                                                                  \
        using Transform = TransformType<double,2>;                                                                          \
        using Point = Kernel<2,double>::Point;                                                                              \
        CIE_TEST_CHECK(SpatialTransform<Transform>);                                                                        \
                                                                                                                            \
        const std::vector<Point> locals {                                                                                   \
            {-1.0, -1.0},                                                                                                   \
            { 1.0,  1.0},                                                                                                   \
                                                                                                                            \
            { 1.0, -1.0},                                                                                                   \
            {-1.0,  1.0},                                                                                                   \
            { 0.0,  0.0}                                                                                                    \
        };                                                                                                                  \
                                                                                                                            \
        const std::vector<Point> transformed {                                                                              \
            { 6.0,  1.0},                                                                                                   \
            { 2.0,  3.0},                                                                                                   \
                                                                                                                            \
            { 2.0,  1.0},                                                                                                   \
            { 6.0,  3.0},                                                                                                   \
            { 4.0,  2.0}                                                                                                    \
        };                                                                                                                  \
                                                                                                                            \
        Transform transform;                                                                                                \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("construction")                                                                              \
            CIE_TEST_CHECK_NOTHROW(transform = Transform(transformed.begin(), transformed.begin() + 2));                    \
        } /*"construction"*/                                                                                                \
                                                                                                                            \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("transformation")                                                                            \
            for (unsigned iPoint=0; iPoint<locals.size(); ++iPoint) {                                                       \
                Point point;                                                                                                \
                transform.evaluate(locals[iPoint], point);                                                                  \
                const auto& rReference = transformed[iPoint];                                                               \
                CIE_TEST_REQUIRE(point.size() == rReference.size());                                                        \
                for (unsigned iComponent=0; iComponent<point.size(); ++iComponent) {                                        \
                    CIE_TEST_CHECK(point[iComponent] == Approx(rReference[iComponent]));                                    \
                }                                                                                                           \
            }                                                                                                               \
        } /*"transformation"*/                                                                                              \
                                                                                                                            \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("derivative")                                                                                \
            decltype(std::declval<Transform>().makeDerivative()) transformDerivative;                                       \
            CIE_TEST_CHECK_NOTHROW(transformDerivative = transform.makeDerivative());                                       \
                                                                                                                            \
            for (unsigned iPoint=0; iPoint<locals.size(); ++iPoint) {                                                       \
                std::array<double,4> jacobian {0.0, 0.0, 0.0, 0.0};                                                         \
                CIE_TEST_REQUIRE(transformDerivative.size() == jacobian.size());                                            \
                transformDerivative.evaluate(locals[iPoint], {jacobian.data(), jacobian.data() + jacobian.size()});         \
                CIE_TEST_CHECK(jacobian[0] == Approx(-2.0));                                                                \
                CIE_TEST_CHECK(jacobian[1] == Approx(0.0).margin(1e-14));                                                   \
                CIE_TEST_CHECK(jacobian[2] == Approx(0.0).margin(1e-14));                                                   \
                CIE_TEST_CHECK(jacobian[3] == Approx(1.0));                                                                 \
            }                                                                                                               \
        } /*"derivative"*/                                                                                                  \
                                                                                                                            \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("inverse")                                                                                   \
            const auto inverseTransform = transform.makeInverse();                                                          \
                                                                                                                            \
            for (unsigned iPoint=0; iPoint<transformed.size(); ++iPoint) {                                                  \
                Point inverse;                                                                                              \
                CIE_TEST_CHECK_NOTHROW(inverseTransform.evaluate(transformed[iPoint], inverse));                            \
                const auto& rReference = locals[iPoint];                                                                    \
                CIE_TEST_REQUIRE(inverse.size() == rReference.size());                                                      \
                for (unsigned iComponent=0; iComponent<inverse.size(); ++iComponent) {                                      \
                    CIE_TEST_CHECK(inverse[iComponent] == Approx(rReference[iComponent]));                                  \
                }                                                                                                           \
            }                                                                                                               \
        } /*"inverse"*/                                                                                                     \
    }


CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST(ScaleTranslateTransform)


CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST(TranslateScaleTransform)


#undef CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST


} // namespace cie::fem::maths
