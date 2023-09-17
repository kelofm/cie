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
        //               |||                                                                                                \
        //               \ /                                                                                                \
        //                v                                                                                                 \
        //                                                                                                                  \
        //   | ( 2, 3)   2---------3   ( 4, 3)                                                                              \
        //   |           |         |                                                                                        \
        //   |           |         |                                                                                        \
        //   |           |         |                                                                                        \
        //   | ( 2, 1)   1---------0   ( 4, 1)                                                                              \
        //   |                                                                                                              \
        // --o------------------------                                                                                      \
        //   |                                                                                                              \
        //   |                                                                                                              \
        */                                                                                                                  \
        using Transform = TransformType<double,2>;                                                                          \
        using Point = Kernel<2,double>::Point;                                                                              \
                                                                                                                            \
        const std::vector<Point> locals {                                                                                   \
            {-1.0, -1.0},                                                                                                   \
            { 1.0,  1.0},                                                                                                   \
                                                                                                                            \
            { 1.0, -1.0},                                                                                                   \
            {-1.0,  1.0}                                                                                                    \
        };                                                                                                                  \
                                                                                                                            \
        const std::vector<Point> transformed {                                                                              \
            { 4.0,  1.0},                                                                                                   \
            { 2.0,  3.0},                                                                                                   \
                                                                                                                            \
            { 2.0,  1.0},                                                                                                   \
            { 4.0,  3.0}                                                                                                    \
        };                                                                                                                  \
                                                                                                                            \
        Transform transform;                                                                                                \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("construction")                                                                              \
            CIE_TEST_CHECK_NOTHROW(transform = Transform(transformed.begin(), transformed.begin() + 1));                    \
        } /*"construction"*/                                                                                                \
                                                                                                                            \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("transformation")                                                                            \
            for (unsigned i_point=0; i_point<locals.size(); ++i_point) {                                                    \
                Point point;                                                                                                \
                transform.evaluate(locals[i_point].data(),                                                                  \
                                   locals[i_point].data() + locals[i_point].size(),                                         \
                                   point.data());                                                                           \
                const auto& r_reference = transformed[i_point];                                                             \
                CIE_TEST_REQUIRE(point.size() == r_reference.size());                                                       \
                for (unsigned i_component=0; i_component<point.size(); ++i_component) {                                     \
                    CIE_TEST_CHECK(point[i_component] == Approx(r_reference[i_component]));                                 \
                }                                                                                                           \
            }                                                                                                               \
        } /*"transformation"*/                                                                                              \
                                                                                                                            \
        {                                                                                                                   \
            CIE_TEST_CASE_INIT("derivative")                                                                                \
            decltype(std::declval<Transform>().makeDerivative()) transformDerivative;                                       \
            CIE_TEST_CHECK_NOTHROW(transformDerivative = transform.makeDerivative());                                       \
                                                                                                                            \
            for (unsigned i_point=0; i_point<locals.size(); ++i_point) {                                                    \
                StaticArray<double,4> jacobian {0.0, 0.0, 0.0, 0.0};                                                        \
                CIE_TEST_REQUIRE(transformDerivative.size() == jacobian.size());                                            \
                transformDerivative.evaluate(locals[i_point].data(),                                                        \
                                             locals[i_point].data() + locals[i_point].size(),                               \
                                             jacobian.data());                                                              \
                CIE_TEST_CHECK(jacobian[0] == Approx(-1.0));                                                                \
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
            for (unsigned i_point=0; i_point<transformed.size(); ++i_point) {                                               \
                Point inverse;                                                                                              \
                CIE_TEST_CHECK_NOTHROW(inverseTransform.evaluate(transformed[i_point].data(),                               \
                                                                 transformed[i_point].data() + transformed[i_point].size(), \
                                                                 inverse.data()));                                          \
                const auto& r_reference = locals[i_point];                                                                  \
                CIE_TEST_REQUIRE(inverse.size() == r_reference.size());                                                     \
                for (unsigned i_component=0; i_component<inverse.size(); ++i_component) {                                   \
                    CIE_TEST_CHECK(inverse[i_component] == Approx(r_reference[i_component]));                               \
                }                                                                                                           \
            }                                                                                                               \
        } /*"inverse"*/                                                                                                     \
    }


CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST(ScaleTranslateTransform)


CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST(TranslateScaleTransform)


#undef CIE_FEM_DEFINE_SCALE_TRANSLATE_TRANSFORM_TEST


} // namespace cie::fem::maths
