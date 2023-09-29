// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/inc/PartitionManager.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"

// --- Utility Includes ---
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/ranges/inc/TransformIterator.hpp"

// --- STL Includes ---
#include <ranges>
#include <iostream>


namespace cie::fem {


CIE_TEST_CASE("PartitionManager", "[graph]")
{
    CIE_TEST_CASE_INIT("PartitionManager")

    AttributeContainer<> empty;
    AttributeContainer<double> root;
    AttributeContainer<int, ParentIndex> mid;
    AttributeContainer<bool, ParentIndex> leaf;

    for (int i=0; i<10; ++i) {
        root.push_back(double(i) / 10.0);
    }
    CIE_TEST_CHECK(root.size() == 10);

    for (std::size_t i=0; i<4; ++i) {
        mid.push_back(10 * i, 2 * i);
    }
    CIE_TEST_CHECK(mid.size() == 4);

    leaf.push_back(false, 3ul);
    leaf.push_back(true, 0ul);
    CIE_TEST_CHECK(leaf.size() == 2);

    // root:   0.0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9
    // mid:    0           10          20          30
    // leaf:   true                                false

    Partition<decltype(root),decltype(empty)> rootPartition(std::move(root), decltype(empty)());
    Partition<decltype(mid),decltype(empty)> midPartition(std::move(mid), decltype(empty)());
    Partition<decltype(leaf),decltype(empty)> leafPartition(std::move(leaf), decltype(empty)());

    {
        const auto [d, i, b] = PartitionManager::collectAttributes<PartitionBase::Item::Vertex>(
            0,
            rootPartition,
            midPartition,
            leafPartition
        );
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, b] = PartitionManager::collectAttributes<PartitionBase::Item::Vertex>(
            1,
            rootPartition,
            midPartition,
            leafPartition
        );
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(b == true);
    }
}


CIE_TEST_CASE("single element", "[graph]")
{
    CIE_TEST_CASE_INIT("single element")

    // Construct bilinear basis functions
    using Basis = maths::AnsatzSpace<maths::Polynomial<Double>,2>;
    auto p_basis = std::make_shared<Basis>(Basis::AnsatzSet {
        maths::Polynomial<Double>({0.5, 0.5}),
        maths::Polynomial<Double>({0.5, -0.5})
    });

    auto p_basisDerivatives = std::make_shared<Basis::Derivative>(p_basis->makeDerivative());

    // Construct mesh
    AttributeContainer<
        unsigned,               // <== vertex ID
        StaticArray<Double,2>   // <== vertex position
    > vertices;

    vertices.push_back(0u, {0.0, 0.0});
    vertices.push_back(1u, {1.0, 1.0});

    AttributeContainer<
        unsigned,                       // <== element ID
        StaticArray<unsigned,2>,        // <== vertex IDs
        decltype(p_basis),              // <== basis functions
        decltype(p_basisDerivatives)    // <== basis function derivatives
    > elements;

    elements.push_back(0u,
                       {0u, 1u},
                       p_basis,
                       p_basisDerivatives);

    //Partition<decltype(vertices),decltype(elements)> partition(std::move(vertices), std::move(elements));

    // Assembly
    Quadrature<Double,2> quadrature(GaussLegendreQuadrature<Double>(2));
    linalg::DynamicEigenMatrix<Double> stiffness(p_basis->size(), p_basis->size());
    std::fill(stiffness.begin(),
              stiffness.end(),
              0.0);

    mp::ThreadLocal<DynamicArray<Double>> tls(DynamicArray<Double>(p_basisDerivatives->size(), 0.0));
    for (Size i_element=0; i_element<elements.size(); ++i_element) {
        elements.visit<StaticArray<unsigned,2>,std::shared_ptr<Basis::Derivative>>(
            i_element,
            [&vertices, &quadrature, &stiffness, &tls](const auto& vertexIDs,
                                                       const auto& p_basisDerivatives){
                const auto vertexCoordinates = std::views::transform(vertexIDs,
                                                                    [&vertices](unsigned vertexID){
                                                                        return vertices.at<StaticArray<Double,2>>(vertexID);
                                                                    });
                const auto transform = maths::ScaleTranslateTransform<Double,2>(
                    vertexCoordinates.begin(),
                    vertexCoordinates.end());
                const Double jacobianDeterminant = transform.makeDerivative().evaluateDeterminant(nullptr, nullptr);
                const auto integrand = maths::makeLambdaExpression<Double>(
                    [&transform, &p_basisDerivatives, &tls, &stiffness](Ptr<const Double> it_argumentBegin,
                                                                        [[maybe_unused]] Ptr<const Double> it_argumentEnd,
                                                                        Ptr<Double> it_output){
                        StaticArray<Double,2> transformedPoint {it_argumentBegin[0], it_argumentBegin[1]};
                        std::cout << transform.size() << std::endl;
                        //transform.evaluate(it_argumentBegin, it_argumentEnd, transformedPoint.data());
                        std::cout << "integration point: " << it_argumentBegin[0] << ", " << it_argumentBegin[1] << "\n"
                                  << "global space:      " << transformedPoint[0] << ", " << transformedPoint[1] << std::endl;
                        auto& r_tls = tls.get<0>();
                        p_basisDerivatives->evaluate(transformedPoint.data(),
                                                     transformedPoint.data() + transformedPoint.size(),
                                                     r_tls.data());
                        Eigen::Map<Eigen::Matrix<Double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> basisDerivativeValues(
                            r_tls.data(),
                            p_basisDerivatives->size() / 2,
                            2);
                        Eigen::Map<Eigen::Matrix<Double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> output(
                            it_output,
                            stiffness.rowSize(),
                            stiffness.columnSize());
                        output = basisDerivativeValues * basisDerivativeValues.transpose();
                    },
                    stiffness.rowSize() * stiffness.columnSize()
                ); // integrand
                quadrature.evaluate(integrand, stiffness.data());
                stiffness.wrapped() *= jacobianDeterminant;
            }
        ); // elements.visit
    } // for i_element

    std::cout << "unconstrained\n" << stiffness.wrapped() << "\n\n";

    for (Size i_constrained : {0ul, 1ul}) {
        for (Size i=0; i<stiffness.rowSize(); ++i) {
            if (i == i_constrained) {
                stiffness(i_constrained, i_constrained) = 1.0;
            } else {
                stiffness(i, i_constrained) = 0.0;
                stiffness(i_constrained, i) = 0.0;
            }
        }
    }
    std::cout << "constrained\n" << stiffness.wrapped() << "\n\n";

    StaticArray<Double,4> displacement {0.0, 0.0, 1.0, 1.0};
    Eigen::Map<Eigen::Matrix<Double,4,1>> disp(displacement.data());
    std::cout << "solution\n" << stiffness.wrapped().fullPivLu().solve(disp) << std::endl;
}


} // namespace cie::fem
