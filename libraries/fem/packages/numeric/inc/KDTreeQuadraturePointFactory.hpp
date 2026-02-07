#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/CompositeDomain.hpp"
#include "packages/numeric/inc/QuadraturePoint.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/utilities/inc/ParametricSpace.hpp"

// --- Geo Includes ---
#include "packages/primitives/inc/Cube.hpp"
#include "packages/trees/inc/ContiguousSpaceTree.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <limits>
#include <span>
#include <algorithm>


namespace cie::fem {


template <
    unsigned PhysicalDimension,
    concepts::Numeric TValue,
    CellLike TCell,
    DomainDataLike TDomainData
>
class KDTreeQuadraturePointFactory {
private:
    static constexpr inline unsigned ParametricDimension = TCell::Dimension;

    using TreePrimitive = geo::Cube<ParametricDimension,TValue>;

    using Tree = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

public:
    using Value = QuadraturePoint<
        ParametricDimension,
        TValue,
        TDomainData>;

    KDTreeQuadraturePointFactory() = default;

    template <CompositeDomainLike<PhysicalDimension,TValue> TDomain>
    KDTreeQuadraturePointFactory(
        Ref<const TDomain> rDomain,
        Ref<const TCell> rCell,
        std::span<const Value> basePoints,
        std::array<unsigned,2> depthRange = {0u, std::numeric_limits<unsigned>::max()})
    {
        CIE_BEGIN_EXCEPTION_TRACING
        TreePrimitive root;
        std::fill_n(
            root.base().data(),
            root.base().size(),
            static_cast<TValue>(-1));
        std::fill_n(
            root.lengths().data(),
            root.lengths().size(),
            static_cast<TValue>(2));
        Tree tree(root);

        const auto visitor = [&, this] (Ref<const typename Tree::Node> rNode, unsigned depth) -> bool {
            // Early exit if the depth is out of range.
            if (depth < depthRange.front()) return true;
            if (depth < depthRange.back()) return false;

            // Fetch the current tree node's geometry (in parametric space).
            std::array<TValue,ParametricDimension> base, lengths;
            tree.getNodeGeometry(
                rNode,
                base.data(),
                lengths.data());

            // Check which subdomain each corner of the current node belong to.
            constexpr unsigned cornerCount = intPow(2, ParametricDimension);
            std::array<TValue,cornerCount*PhysicalDimension> physicalCorners;
            std::array<typename TDomain::DomainData,cornerCount> subdomains;

            Ptr<TValue> pCornerBegin = physicalCorners.data();
            ParametricSpace<ParametricDimension,TValue,ParametricSpaceType::Cartesian>::iterateCorners(
                [&subdomains, &pCornerBegin, &rNode, &rDomain, &rCell, &base, &lengths] (std::span<const std::uint8_t,cornerCount> state) {
                    // Compute the tree node's corner (in parametric space).
                    std::array<TValue,cornerCount> parametricCorner;
                    for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
                        parametricCorner[iDimension] = base[iDimension];
                        if (state[iDimension]) parametricCorner[iDimension] += lengths[iDimension];
                    } // for iDimension in range(ParametricDimension)

                    // Transform the corner from parametric to physical space.
                    rCell.spatialTransform().evaluate(
                        parametricCorner,
                        std::span<TValue>(pCornerBegin, PhysicalDimension));

                    pCornerBegin += PhysicalDimension;
                } // iterateCorners lambda
            ); // ParametricSpace::iterateCorners

            // Decide whether to generate quadrature points.
            // - If all corners lie in the same domain, quadrature points can
            //   be immediately generated and the branch need not be
            //   further explored.
            // - If not all corners belong to the same domain, further
            //   partitioning is necessary.
            rDomain.whichSubdomain(
                std::span<const TValue>(physicalCorners),
                std::span<typename TDomain::DomainData>(subdomains));
            const bool isHomogeneousSubdomain = std::adjacent_find(
                subdomains.begin(),
                subdomains.end(),
                std::not_equal_to<typename TDomain::DomainData>()
            ) == subdomains.end();

            if (isHomogeneousSubdomain && subdomains.front()) {
                // Compute the spatial transform that maps from the node's
                // local space to the cell's parametric space.
                std::array<
                    std::array<
                        TValue,
                        ParametricDimension>
                    ,2
                > transformedCorners;
                std::copy_n(
                    base.data(),
                    ParametricDimension,
                    transformedCorners.front().data());
                std::transform(
                    base.begin(),
                    base.end(),
                    lengths.begin(),
                    transformedCorners.back().begin(),
                    std::plus<TValue>());
                maths::ScaleTranslateTransform<TValue,ParametricDimension> transform(
                    transformedCorners.data(),
                    transformedCorners.data() + 1);
                const TValue inverseJacobianDeterminant = transform.makeDerivative().evaluateDeterminant(std::span<TValue>());

                _points.reserve(_points.size() + basePoints.size());
                for (Ref<const Value> rBasePoint : basePoints) {
                    Ref<Value> rPoint = _points.emplace_back(rBasePoint);
                    transform.evaluate(
                        rBasePoint.position(),
                        rPoint.position());
                    rPoint.weight() *= inverseJacobianDeterminant;
                }  // for rBasePoint in basePoints
            } // if isHomogeneousSubdomain

            return !isHomogeneousSubdomain;
        }; // visitor

        tree.scan(visitor);
        CIE_END_EXCEPTION_TRACING
    }

    unsigned operator()(std::span<Value> out) {
        const std::size_t copyCount = std::min<std::size_t>(
            std::distance(_pBegin, _points.data() + _points.size()),
            out.size());
        std::copy_n(
            _points.data(),
            copyCount,
            out.data());
        _pBegin += copyCount;
        return copyCount;
    }

private:
    DynamicArray<Value> _points;

    Ptr<const Value> _pBegin;
}; // class KDTreeQuadraturePointFactory


} // namespace cie::fem
