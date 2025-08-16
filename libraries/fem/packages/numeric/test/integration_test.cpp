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

    const auto domain = [] (std::span<const double> in) -> bool {
        return linalg::norm2(in) < 1;
    };

    const auto treePredicate = [&tree, &domain] (Ref<const QuadTree::Node> rNode, unsigned level) -> bool {
        if (10 < level) {return false;}
        QuadTree::Point base;
        double edge;
        tree.getNodeGeometry(rNode, base.data(), &edge);
        bool isInside = domain(base);
        base[0] += edge;
        if (isInside != domain(base)) return true;
        base[1] += edge;
        if (isInside != domain(base)) return true;
        base[0] -= edge;
        if (isInside != domain(base)) return true;
        return false;
    };
    tree.scan(treePredicate);

    // Construct an integrand that will be evaluated over the domain
    const auto integrand = maths::makeLambdaExpression<double>([&domain] (std::span<const double> in,
                                                                          std::span<double> out) {
        if (domain(in)) {
            // A unit halfsphere
            out.front() = std::sqrt(1.0 - std::pow(in[0], 2) - std::pow(in[1], 2));
        } else {
            out.front() = 0;
        }
    }, 1);

    double integral = 0;
    for (const auto& rNode : tree) {
        if (rNode.isLeaf()) {
            // Recover the node's geometry
            double base[2];
            double edge;
            tree.getNodeGeometry(rNode, base, &edge);

            // Construct the transformation that maps
            // local space [-1,1]^2 to the node's geometry
            StaticArray<double,2> transformed[2];
            transformed[0] = {base[0], base[1]};
            transformed[1] = {base[0] + edge, base[1] + edge};
            const maths::ScaleTranslateTransform<double,2> transform(transformed, transformed + 2);
            const double determinant = transform.makeDerivative().evaluateDeterminant({&edge, 0});

            // Construct the transformed integrand
            const auto transformedIntegrand = maths::makeLambdaExpression<double>(
                [&transform, determinant, &integrand] (std::span<const double> in,
                                                       std::span<double> out) {
                StaticArray<double,2> transformedPoint;
                transform.evaluate(in, transformedPoint);
                integrand.evaluate(transformedPoint, out);
                out.front() *= determinant;
            }, integrand.size());

            double term;
            quadrature.evaluate(transformedIntegrand, {&term, 1});
            integral += term;
        } // if rNode.isLeaf()
    } // for node in tree

    CIE_TEST_CHECK(integral == Approx(/*sphere volume*/ 4.0 / 3.0 * std::numbers::pi /*but only an eighth*/ / 8));
}


} // namespace cie::fem
