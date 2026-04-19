#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/CompositeDomain.hpp"
#include "packages/numeric/inc/QuadraturePoint.hpp"
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
    DomainDataLike TDomainData = void>
class KDTreeQuadraturePointFactory {
private:
    static constexpr inline unsigned ParametricDimension = TCell::Dimension;

    using TreePrimitive = geo::Cube<ParametricDimension,TValue>;

    using Tree = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

public:
    static constexpr inline unsigned Dimension = ParametricDimension;

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
        std::array<unsigned,2> depthRange = {0u, std::numeric_limits<unsigned>::max()}) {
            CIE_BEGIN_EXCEPTION_TRACING
            TreePrimitive root;
            std::fill_n(
                root.base().data(),
                root.base().size(),
                static_cast<TValue>(-1));
            root.length() = static_cast<TValue>(2);
            Tree tree(root);

            std::vector<TValue> buffer;
            std::vector<TDomainData> subdomains;
            std::vector<TValue> physicalCoordinates;

            const auto visitor = [&, this] (Ref<const typename Tree::Node> rNode, unsigned depth) -> bool {
                // Early exit if the depth is out of range.
                if (depth < depthRange.front()) return true;
                if (depthRange.back() < depth) return false;

                // Fetch the current tree node's geometry (in parametric space).
                std::array<TValue,ParametricDimension> base;
                TValue length = 0;
                tree.getNodeGeometry(
                    rNode,
                    base.data(),
                    &length);

                //if ((isHomogeneousSubdomain && subdomains.front()) || depth == depthRange.back()) {
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
                    transformedCorners.back().begin(),
                    [length] (TValue component) -> TValue {return component + length;});
                const maths::ScaleTranslateTransform<TValue,ParametricDimension> localTransform(
                    transformedCorners.data(),
                    transformedCorners.data() + 2);
                const auto localJacobian = localTransform.makeDerivative();

                // Transform quadrature point locations from the standard
                // quadrature domain to the cell's parametric domain.
                // At the same time, compute the physical location too
                // to later determine which subdomain each integration point
                // lies in.
                _points.reserve(_points.size() + basePoints.size());
                physicalCoordinates.resize(basePoints.size() * ParametricDimension);
                buffer.resize(std::max(
                    localTransform.bufferSize(),
                    rCell.makeSpatialTransform().bufferSize()));

                for (std::size_t iPoint=0ul; iPoint<basePoints.size(); ++iPoint) {
                    // Make a copy of the integration rule's template.
                    auto point = basePoints[iPoint];

                    // Map the quadrature point to the cell's parametric space.
                    localTransform.evaluate(
                        Kernel<ParametricDimension,TValue>::decay(basePoints[iPoint].position()),
                        Kernel<ParametricDimension,TValue>::decay(point.position()),
                        buffer);

                    // Map the quadrature point to physical space.
                    std::span<TValue,TCell::PhysicalDimension> physicalPosition(
                        physicalCoordinates.data() + iPoint * TCell::PhysicalDimension,
                        TCell::PhysicalDimension);
                    rCell.transform(
                        point.position(),
                        Kernel<TCell::PhysicalDimension,TValue>::template cast<PhysicalCoordinate<TValue>>(physicalPosition),
                        buffer);

                    // Compute the jacobian's determinant and scale the quadrature
                    // point's weight with it.
                    const TValue jacobianDeterminant = localJacobian.evaluateDeterminant(
                        Kernel<ParametricDimension,TValue>::decay(point.position()),
                        buffer);
                    point.weight() *= jacobianDeterminant;

                    // Register the quadrature point.
                    _points.emplace_back(point);
                }  // for rBasePoint in basePoints

                // Find subdomains.
                subdomains.resize(basePoints.size());
                rDomain.subdomain(
                    std::span<const TValue>(physicalCoordinates),
                    std::span<typename TDomain::DomainData>(subdomains));

                const bool isHomogeneousSubdomain = std::adjacent_find(
                    subdomains.begin(),
                    subdomains.end(),
                    std::not_equal_to<typename TDomain::DomainData>()
                ) == subdomains.end();

                // Decide whether to keep quadrature points.
                bool keep = true;
                keep &= (isHomogeneousSubdomain && subdomains.front());
                keep |= depth == depthRange.back();
                keep |= depth == depthRange.front() && isHomogeneousSubdomain;

                if (keep) {
                    // Keep the generated quadrature points and set their subdomains.
                    for (std::size_t iPoint=0ul; iPoint<basePoints.size(); ++iPoint) {
                        _points[_points.size() - basePoints.size() + iPoint].data() = subdomains[iPoint];
                    }
                } /*if isHomogeneousSubdomain*/ else {
                    // Discard generated quadrature points.
                    _points.resize(_points.size() - basePoints.size());
                }

                return !isHomogeneousSubdomain;
            }; // visitor

            tree.scan(visitor);
            _pBegin = _points.data();
            CIE_END_EXCEPTION_TRACING
    }

    unsigned operator()(std::span<Value> out) {
        const std::size_t copyCount = std::min<std::size_t>(
            std::distance<Ptr<const Value>>(_pBegin, _points.data() + _points.size()),
            out.size());
        std::copy_n(
            _pBegin,
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
