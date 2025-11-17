// --- External Includes ---
#include "pybind11/pybind11.h"


// --- FEM Includes ---
#include "packages/maths/inc/OrthogonalScaleTransform.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/maths/inc/AffineTransform.hpp"
#include "packages/maths/inc/ProjectiveTransform.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "python/bindings/inc/stl_extension.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <numeric>


namespace cie::fem::maths {


/** @brief Add python interface for vector => vector @ref Expression s.
 *  @tparam TExpression Class to add member functions from, satisfying @ref Expression.
 *  @tparam InDimension Input dimension.
 *  @tparam OutDimension Output dimension.
 *  @tparam TPyBindClass Pybind11 wrapper for @a TExpression.
 *  @param rClass Reference to the pybind wrapper to add the expression interface to.
 */
template <Expression TExpression,
          unsigned InDimension,
          unsigned OutDimension,
          class TPyBindClass>
void addVectorExpressionInterface(TPyBindClass rClass)
{
    using TValue = typename TExpression::Value;
    rClass.def("size",
                &TExpression::size,
                "Get the number of scalar components returned by the 'evaluate' member function.")
    .def("evaluate", [] (Ref<const TExpression> rThis,
                         Ref<const StaticArray<TValue,InDimension>> rArguments) {
        StaticArray<Size,1> array;
        array.front() = OutDimension;
        auto adaptor = makeNumpyArray<TValue,1>(array.data());
        rThis.evaluate(rArguments, {adaptor.mutable_data(), static_cast<std::size_t>(adaptor.size())});
        return adaptor;
    });
}


/** @brief Add python interface for vector => vector @ref Expression s.
 *  @tparam TExpression Class to add member functions from, satisfying @ref Expression.
 *  @tparam InDimension Input dimension.
 *  @tparam OutRows Number of rows in the output.
 *  @tparam OutColumns Number of columns in the output.
 *  @tparam TPyBindClass Pybind11 wrapper for @a TExpression.
 *  @param rClass Reference to the pybind wrapper to add the expression interface to.
 */
template <Expression TExpression,
          Size InDimension,
          Size OutRows,
          Size OutColumns,
          class TPyBindClass>
void addMatrixExpressionInterface(TPyBindClass rClass)
{
    using TValue = typename TExpression::Value;
    rClass.def("size",
                &TExpression::size,
                "Get the number of scalar components returned by the 'evaluate' member function.")
    .def("evaluate", [] (Ref<const TExpression> rThis,
                         Ref<const StaticArray<TValue,InDimension>> rArguments) {
        auto array = makeNumpyArray<TValue,2>(StaticArray<Size,2> {OutRows, OutColumns} .data());
        rThis.evaluate(rArguments, {array.mutable_data(), static_cast<std::size_t>(array.size())});
        return array;
    });
}


void makeFEMMathsBindings(Ref<pybind11::module_> rModule)
{
    auto submodule = rModule.def_submodule("maths", "");

    {
        auto orthogonalScaleTransformDerivative2D = pybind11::class_<OrthogonalScaleTransformDerivative<double,2>>(
            submodule,
            "OrthogonalScaleTransformDerivative2D"
        )   .def(pybind11::init<>())
            .def("evaluateDeterminant", [] (Ref<OrthogonalScaleTransformDerivative<double,2>> rThis, Ref<const StaticArray<double,2>> rArgument) {
                return rThis.evaluateDeterminant(rArgument);
            })
            ;
        addMatrixExpressionInterface<OrthogonalScaleTransformDerivative<double,2>,
                                    /*InDimension*/ 2,
                                    /*OutRows    */ 2,
                                    /*OutColumns */ 2>(orthogonalScaleTransformDerivative2D);
    }

    {
        auto orthogonalScaleTransform2D = pybind11::class_<OrthogonalScaleTransform<double,2>>(
            submodule,
            "OrthogonalScaleTransform2D"
        )   .def(pybind11::init<>())
            .def(pybind11::init([] (Ref<const DynamicArray<StaticArray<double,2>>> rTransformed) {
                return OrthogonalScaleTransform<double,2>(rTransformed.begin(), rTransformed.end());
            }))
            .def("makeDerivative", &OrthogonalScaleTransform<double,2>::makeDerivative)
            .def("makeInverse", &OrthogonalScaleTransform<double,2>::makeInverse)
            ;
        addVectorExpressionInterface<OrthogonalScaleTransform<double,2>,
                                    /*InDimension */ 2,
                                    /*OutDimension*/ 2>(orthogonalScaleTransform2D);
    }

    {
        auto scaleTranslateTransformDerivative2D = pybind11::class_<ScaleTranslateTransformDerivative<double,2>>(
            submodule,
            "ScaleTranslateTransformDerivative2D"
        )   .def(pybind11::init<>())
            .def("evaluateDeterminant", [] (Ref<ScaleTranslateTransformDerivative<double,2>> rThis, Ref<const StaticArray<double,2>> rArgument) {
                return rThis.evaluateDeterminant(rArgument);
            })
            ;
        addMatrixExpressionInterface<ScaleTranslateTransformDerivative<double,2>,
                                     2,
                                     2,
                                     2>(scaleTranslateTransformDerivative2D);
    }

    {
        auto scaleTranslateTransform2D = pybind11::class_<ScaleTranslateTransform<double,2>>(
            submodule,
            "ScaleTranslateTransform2D"
        )   .def(pybind11::init<>())
            .def(pybind11::init([] (Ref<const DynamicArray<StaticArray<double,2>>> rTransformed) {
                return ScaleTranslateTransform<double,2>(rTransformed.begin(), rTransformed.end());
            }))
            .def("makeDerivative", &ScaleTranslateTransform<double,2>::makeDerivative)
            .def("makeInverse", &ScaleTranslateTransform<double,2>::makeInverse)
            ;
        addVectorExpressionInterface<ScaleTranslateTransform<double,2>,
                                     2,
                                     2>(scaleTranslateTransform2D);
    }

    {
        auto translateScaleTransform2D = pybind11::class_<TranslateScaleTransform<double,2>>(
            submodule,
            "TranslateScaleTransform2D"
        )   .def(pybind11::init<>())
            .def(pybind11::init([] (Ref<const DynamicArray<StaticArray<double,2>>> rTransformed) {
                return TranslateScaleTransform<double,2>(rTransformed.begin(), rTransformed.end());
            }))
            .def("makeDerivative", &TranslateScaleTransform<double,2>::makeDerivative)
            .def("makeInverse", &TranslateScaleTransform<double,2>::makeInverse)
            ;
        addVectorExpressionInterface<TranslateScaleTransform<double,2>,
                                     2,
                                     2>(translateScaleTransform2D);
    }

    {
        auto affineTransformDerivative2D = pybind11::class_<AffineTransformDerivative<double,2>>(submodule, "AffineTransformDerivative2D")
            .def(pybind11::init<>())
            .def("evaluateDeterminant", [] (Ref<AffineTransformDerivative<double,2>> rThis, Ref<const StaticArray<double,2>> rArgument) {
                return rThis.evaluateDeterminant(rArgument);
            })
            ;
        addMatrixExpressionInterface<AffineTransformDerivative<double,2>,
                                    /*InDimension*/ 2,
                                    /*OutRows    */ 2,
                                    /*OutColumns */ 2>(affineTransformDerivative2D);
    }

    {
        auto affineTransform2D = pybind11::class_<AffineTransform<double,2>>(submodule, "AffineTransform2D")
            .def(pybind11::init<>())
            .def(pybind11::init([](Ref<const DynamicArray<StaticArray<double,2>>> rTransformed) {
                    StaticArray<AffineTransform<double,2>::Point,3> transformed;
                    CIE_OUT_OF_RANGE_CHECK(rTransformed.size() == transformed.size());

                    for (unsigned iPoint=0; iPoint<transformed.size(); ++iPoint) {
                        std::copy(rTransformed[iPoint].begin(),
                                  rTransformed[iPoint].end(),
                                  transformed[iPoint].begin());
                    }

                    return AffineTransform<double,2>(transformed);
                }))
            .def("makeDerivative", &AffineTransform<double,2>::makeDerivative)
            .def("makeInverse", &AffineTransform<double,2>::makeInverse)
            ;
        addVectorExpressionInterface<AffineTransform<double,2>,
                                    /*InDimension */ 2,
                                    /*OutDimension*/ 2>(affineTransform2D);
    }

    {
        auto projectiveTransformDerivative2D = pybind11::class_<ProjectiveTransformDerivative<double,2>>(submodule, "ProjectiveTransformDerivative2D")
            .def(pybind11::init<>())
            .def("evaluateDeterminant", [] (Ref<ProjectiveTransformDerivative<double,2>> rThis, Ref<const StaticArray<double,2>> rArgument) {
                return rThis.evaluateDeterminant(rArgument);
            })
            ;
        addMatrixExpressionInterface<ProjectiveTransformDerivative<double,2>,
                                    /*InDimension*/ 2,
                                    /*OutRows    */ 2,
                                    /*OutColumns */ 2>(projectiveTransformDerivative2D);
    }

    {
        auto projectiveTransform2D = pybind11::class_<ProjectiveTransform<double,2>>(submodule, "ProjectiveTransform2D")
            .def(pybind11::init<>())
            .def(pybind11::init([](const std::vector<StaticArray<double,2>>& rTransformed){
                    std::vector<ProjectiveTransform<double,2>::Point> transformed;
                    transformed.reserve(rTransformed.size());
                    for (const auto& rVertex : rTransformed) {
                        transformed.push_back({rVertex[0], rVertex[1]});
                    }
                    return ProjectiveTransform<double,2>(transformed);
                }))
            .def("makeDerivative", &ProjectiveTransform<double,2>::makeDerivative)
            .def("makeInverse", &ProjectiveTransform<double,2>::makeInverse)
            ;
        addVectorExpressionInterface<ProjectiveTransform<double,2>,
                                    /*InDimension */ 2,
                                    /*OutDimension*/ 2>(projectiveTransform2D);
    }
}


} // namespace cie::fem::maths
