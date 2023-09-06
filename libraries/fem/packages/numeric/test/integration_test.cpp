// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- GEO Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/primitives/inc/Cube.hpp"
#include "packages/trees/inc/ContiguousSpaceTree.hpp"

// --- Linalg Includes ---
#include "packages/types/impl/typeoperations_impl.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- STL Includes ---
#include <numbers>
#include <cmath>


namespace cie::fem {


CIE_TEST_CASE("integration", "[fem]")
{
    CIE_TEST_CASE_INIT("integration")

    using QuadTree = geo::ContiguousSpaceTree<geo::Cube<2,double>,unsigned>;

    const Quadrature<double,2> quadrature(GaussLegendreQuadrature<double>(8));
    QuadTree tree(QuadTree::Point {0.0, 0.0}, 1.0);

    const auto domain = [] (Ptr<const double> it_begin,
                            Ptr<const double>) -> bool {
        StaticArray<double,2> point {*it_begin, *(it_begin + 1)};
        return linalg::norm2(point) < 1;
    };

    const auto treePredicate = [&tree, &domain] (Ref<const QuadTree::Node> r_node) -> bool {
        QuadTree::Point base;
        double edge;
        tree.getNodeGeometry(r_node, base.data(), &edge);
        bool isInside = domain(base.begin(), base.end());
        base[0] += edge;
        if (isInside != domain(base.begin(), base.end())) return true;
        base[1] += edge;
        if (isInside != domain(base.begin(), base.end())) return true;
        base[0] -= edge;
        if (isInside != domain(base.begin(), base.end())) return true;
        return false;
    };
    tree.scan(treePredicate, 10);

    // Construct an integrand that will be evaluated over the domain
    const auto integrand = maths::makeLambdaExpression<double>([&domain] (Ptr<const double> it_begin,
                                                                          Ptr<const double> it_end,
                                                                          Ptr<double> it_out) {
        if (domain(it_begin, it_end)) {
            // A unit halfsphere
            *it_out = std::sqrt(1.0 - std::pow(*it_begin, 2) - std::pow(*(it_begin + 1), 2));
        } else {
            *it_out = 0;
        }
    }, 1);

    double integral = 0;
    for (const auto& r_node : tree) {
        if (r_node.isLeaf()) {
            // Recover the node's geometry
            double base[2];
            double edge;
            tree.getNodeGeometry(r_node,
                                 base,
                                 &edge);

            // Construct the transformation that maps
            // local space [-1,1]^2 to the node's geometry
            StaticArray<double,2> transformed[2];
            transformed[0] = {base[0], base[1]};
            transformed[1] = {base[0] + edge, base[1] + edge};
            const maths::ScaleTranslateTransform<double,2> transform(transformed, transformed + 2);
            const double determinant = transform.makeDerivative().evaluateDeterminant(nullptr, nullptr);

            // Construct the transformed integrand
            const auto transformedIntegrand = maths::makeLambdaExpression<double>(
                [&transform, determinant, &integrand] (Ptr<const double> it_begin,
                                                       Ptr<const double> it_end,
                                                       Ptr<double> it_out) {
                StaticArray<double,2> transformedPoint;
                transform.evaluate(it_begin, it_end, transformedPoint.data());
                integrand.evaluate(transformedPoint.begin(),
                                   transformedPoint.end(),
                                   it_out);
                *it_out *= determinant;
            }, integrand.size());

            double term;
            quadrature.evaluate(transformedIntegrand, &term);
            integral += term;
        } // if r_node.isLeaf()
    } // for node in tree

    CIE_TEST_CHECK(integral == Approx(/*sphere volume*/ 4.0 / 3.0 * std::numbers::pi /*but only an eight*/ / 8));
}


} // namespace cie::fem
