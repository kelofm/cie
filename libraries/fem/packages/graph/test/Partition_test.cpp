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
    auto pBasis = std::make_shared<Basis>(Basis::AnsatzSet {
        maths::Polynomial<Double>({0.5, 0.5}),
        maths::Polynomial<Double>({0.5, -0.5})
    });

    auto pBasisDerivatives = std::make_shared<Basis::Derivative>(pBasis->makeDerivative());

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
        decltype(pBasis),              // <== basis functions
        decltype(pBasisDerivatives)    // <== basis function derivatives
    > elements;

    elements.push_back(0u,
                       {0u, 1u},
                       pBasis,
                       pBasisDerivatives);

    //Partition<decltype(vertices),decltype(elements)> partition(std::move(vertices), std::move(elements));

    // Assembly
    Quadrature<Double,2> quadrature(GaussLegendreQuadrature<Double>(2));
    linalg::DynamicEigenMatrix<Double> stiffness(pBasis->size(), pBasis->size());
    std::fill(stiffness.begin(),
              stiffness.end(),
              0.0);

    mp::ThreadLocal<DynamicArray<Double>> tls(DynamicArray<Double>(pBasisDerivatives->size(), 0.0));
    for (Size iElement=0; iElement<elements.size(); ++iElement) {
        elements.visit<StaticArray<unsigned,2>,std::shared_ptr<Basis::Derivative>>(
            iElement,
            [&vertices, &quadrature, &stiffness, &tls](const auto& vertexIDs,
                                                       const auto& pBasisDerivatives){
                const auto vertexCoordinates = std::views::transform(vertexIDs,
                                                                    [&vertices](unsigned vertexID){
                                                                        return vertices.at<StaticArray<Double,2>>(vertexID);
                                                                    });
                const auto transform = maths::ScaleTranslateTransform<Double,2>(
                    vertexCoordinates.begin(),
                    vertexCoordinates.end());
                const Double jacobianDeterminant = transform.makeDerivative().evaluateDeterminant(nullptr, nullptr);
                const auto integrand = maths::makeLambdaExpression<Double>(
                    [&transform, &pBasisDerivatives, &tls, &stiffness](Ptr<const Double> itArgumentBegin,
                                                                        [[maybe_unused]] Ptr<const Double> itArgumentEnd,
                                                                        Ptr<Double> itOutput){
                        StaticArray<Double,2> transformedPoint {itArgumentBegin[0], itArgumentBegin[1]};
                        std::cout << transform.size() << std::endl;
                        //transform.evaluate(itArgumentBegin, itArgumentEnd, transformedPoint.data());
                        std::cout << "integration point: " << itArgumentBegin[0] << ", " << itArgumentBegin[1] << "\n"
                                  << "global space:      " << transformedPoint[0] << ", " << transformedPoint[1] << std::endl;
                        auto& rTls = tls.get<0>();
                        pBasisDerivatives->evaluate(transformedPoint.data(),
                                                     transformedPoint.data() + transformedPoint.size(),
                                                     rTls.data());
                        Eigen::Map<Eigen::Matrix<Double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> basisDerivativeValues(
                            rTls.data(),
                            pBasisDerivatives->size() / 2,
                            2);
                        Eigen::Map<Eigen::Matrix<Double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> output(
                            itOutput,
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
    } // for iElement

    std::cout << "unconstrained\n" << stiffness.wrapped() << "\n\n";

    for (Size iConstrained : {0ul, 1ul}) {
        for (Size i=0; i<stiffness.rowSize(); ++i) {
            if (i == iConstrained) {
                stiffness(iConstrained, iConstrained) = 1.0;
            } else {
                stiffness(i, iConstrained) = 0.0;
                stiffness(iConstrained, i) = 0.0;
            }
        }
    }
    std::cout << "constrained\n" << stiffness.wrapped() << "\n\n";

    StaticArray<Double,4> displacement {0.0, 0.0, 1.0, 1.0};
    Eigen::Map<Eigen::Matrix<Double,4,1>> disp(displacement.data());
    std::cout << "solution\n" << stiffness.wrapped().fullPivLu().solve(disp) << std::endl;
}


} // namespace cie::fem
