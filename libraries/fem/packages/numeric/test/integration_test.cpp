// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/Serializer.hpp"

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
#include "packages/utilities/inc/kernel.hpp"

// --- STL Includes ---
#include <numbers>
#include <cmath>


namespace cie::fem {


bool testDomain(std::span<const double> in) {
    return linalg::norm2(in) < 1;
}


struct IntegrationTestIntegrand
    :   public maths::ExpressionTraits<double>,
        public cie::io::TriviallySerializableBase {
    IntegrationTestIntegrand() {}

    IntegrationTestIntegrand(Ref<const maths::ScaleTranslateTransform<double,2>> rTransform)
        : _transform(rTransform) {}

    void evaluate(
        typename maths::ExpressionTraits<double>::ConstSpan in,
        typename maths::ExpressionTraits<double>::Span out,
        typename maths::ExpressionTraits<double>::BufferSpan buffer) const {
            CIE_TEST_REQUIRE(in.size() == 2);
            CIE_TEST_REQUIRE(out.size() == this->size());
            CIE_TEST_REQUIRE(this->bufferSize() <= buffer.size());
            std::array<double,2> transformed;
            _transform.evaluate(in, transformed, buffer);
            if (testDomain(transformed)) {
                // A unit halfsphere
                out.front() = std::sqrt(1.0 - std::pow(transformed[0], 2) - std::pow(transformed[1], 2));
                out.front() *= _transform.makeDerivative().evaluateDeterminant(in, buffer);
            } else {
                out.front() = 0;
            }
    }

    unsigned size() const {
        return 1u;
    }

    unsigned bufferSize() const {
        return std::max<unsigned>(
            _transform.bufferSize(),
            _transform.makeDerivative().bufferSize());
    }

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary) const {
            cie::io::BinarySerializer::serialize(
                rStream,
                _transform);
    }

    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<IntegrationTestIntegrand> rInstance,
        tags::Binary) {
            cie::io::BinarySerializer::deserialize(
                rStream,
                rInstance._transform);
    }

    maths::ScaleTranslateTransform<double,2> _transform;
}; // struct IntegrationTestIntegrand


CIE_TEST_CASE("integration", "[fem]") {
    CIE_TEST_CASE_INIT("integration")

    using QuadTree = geo::ContiguousSpaceTree<geo::Cube<2,double>,unsigned>;

    const Quadrature<double,2> quadrature(GaussLegendreQuadrature<double>(8));
    QuadTree tree(QuadTree::Point {0.0, 0.0}, 1.0);

    const auto treePredicate = [&tree] (Ref<const QuadTree::Node> rNode, unsigned level) -> bool {
        if (10 < level) {return false;}
        QuadTree::Point base;
        double edge;
        tree.getNodeGeometry(rNode, base.data(), &edge);
        bool isInside = testDomain(base);
        base[0] += edge;
        if (isInside != testDomain(base)) return true;
        base[1] += edge;
        if (isInside != testDomain(base)) return true;
        base[0] -= edge;
        if (isInside != testDomain(base)) return true;
        return false;
    };
    tree.scan(treePredicate);

    double integral = 0;
    std::vector<double> buffer;

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
            const maths::ScaleTranslateTransform<double,2> transform(
                transformed,
                transformed + 2);

            // Construct the transformed integrand
            const IntegrationTestIntegrand integrand(transform);

            double term;
            buffer.resize(quadrature.bufferSize(integrand));
            quadrature.evaluate(
                integrand,
                buffer,
                {&term, 1});
            integral += term;
        } // if rNode.isLeaf()
    } // for node in tree

    CIE_TEST_CHECK(integral == Approx(/*sphere volume*/ 4.0 / 3.0 * std::numbers::pi /*but only an eighth*/ / 8));
}


} // namespace cie::fem
