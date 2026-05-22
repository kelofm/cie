// --- Internal Includes ---
#include "embeddedPoisson2D/MeshData.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"


namespace cie::fem {


MeshData::MeshData()
    : MeshBase<Ansatz>(),
      _quadraturePointSet()
{}


MeshData::MeshData(
    RightRef<Ansatz> rAnsatzSpace,
    RightRef<std::vector<std::pair<DomainData,std::vector<Scalar>>>> rDomainTriangles,
    std::span<const std::pair<DomainData,Scalar>> domainMap,
    Ref<const cie::io::JSONObject> rConfiguration)
        :   MeshBase<Ansatz>(std::span<const Ansatz>(&rAnsatzSpace, 1)),
            _quadraturePointSet(),
            _domainTriangles(std::move(rDomainTriangles)),
            _domainMap(domainMap.begin(), domainMap.end()),
            _cells(),
            _treeDepthRange() {
                // Cache quadrature points.
                // Generate 1D quadrature points.
                DynamicArray<QuadraturePoint<1,Scalar,Scalar>> basePoints;
                OuterProductQuadraturePointFactory<Dimension,Scalar,Scalar> generator;

                CIE_BEGIN_EXCEPTION_TRACING
                    const std::size_t integrationOrder = rConfiguration["integration"]["template"]["order"].as<std::size_t>();
                    _treeDepthRange.front() = rConfiguration["integration"]["min-tree-depth"].as<std::size_t>();
                    _treeDepthRange.back() = rConfiguration["integration"]["max-tree-depth"].as<std::size_t>();

                    GaussLegendreQuadrature<Scalar> quadrature(integrationOrder);
                    basePoints.reserve(quadrature.numberOfNodes());
                    for (std::size_t iNode=0ul; iNode<quadrature.numberOfNodes(); ++iNode) {
                        basePoints.emplace_back(
                            quadrature.nodes()[iNode],
                            quadrature.weights()[iNode]);
                    }
                    generator = OuterProductQuadraturePointFactory<Dimension,Scalar,Scalar>(basePoints);
                CIE_END_EXCEPTION_TRACING

                // Generate nD quadrature points.
                constexpr std::size_t initialSize = 0x10;
                constexpr Scalar growFactor = 2.0;
                std::size_t pointCount = 0ul;
                _quadraturePointSet.resize(initialSize);

                while (true) {
                    // Extend the container if necessary.
                    if (!(pointCount < _quadraturePointSet.size())) {
                        _quadraturePointSet.resize(growFactor * _quadraturePointSet.size());
                    }

                    std::span<QuadraturePoint<Dimension,Scalar,Scalar>> targetSpan(
                        _quadraturePointSet.data() + pointCount,
                        _quadraturePointSet.data() + _quadraturePointSet.size());

                    // Request a new batch of quadrature points.
                    const unsigned newPointCount = generator(targetSpan);
                    pointCount += newPointCount;

                    if (!newPointCount) break;
                }

                _quadraturePointSet.resize(pointCount);
                _quadraturePointSet.shrink_to_fit();
}


KDTreeQuadraturePointFactory<
    Dimension,
    Scalar,
    Cell,
    Scalar
> MeshData::makeQuadratureRule(Ref<const Cell> rCell) const {
    return KDTreeQuadraturePointFactory<
        Dimension,
        Scalar,
        Cell,
        Scalar>(
            *this,
            rCell,
            std::span<const QuadraturePoint<Dimension,Scalar,Scalar>>(
                _quadraturePointSet.data(),
                _quadraturePointSet.size()),
        _treeDepthRange
    );
}


namespace impl {


/// @details Checks whether a reference point is on the negative side
///          (opposite the normal) of a line segment defined by two points.
constexpr bool negativeSide(
    std::span<const Scalar,Dimension> reference,
    std::span<const Scalar,Dimension> begin,
    std::span<const Scalar,Dimension> end) noexcept {
        std::array<Scalar,Dimension> normal, relative;

        std::transform(
            reference.begin(),
            reference.end(),
            begin.begin(),
            relative.begin(),
            std::minus<Scalar>());

        if constexpr (Dimension == 2) {
            const std::array<Scalar,2> diff {
                end[0] - begin[0],
                end[1] - begin[1]};
            normal[0] = diff[1];
            normal[1] = -diff[0];
        } else static_assert(Dimension == 2, "unsupported dimension");

        return std::inner_product(
            relative.begin(),
            relative.end(),
            normal.begin(),
            static_cast<Scalar>(0)) < static_cast<Scalar>(0);
}


constexpr bool isInTriangle(
    std::span<const Scalar,Dimension> reference,
    std::span<const Scalar,3*Dimension> triangle) noexcept {
        bool output = true;
        for (std::size_t iSegment=0ul; iSegment<3ul; ++iSegment) {
            std::array<Scalar,2*Dimension> segment;
            for (std::size_t iVertex=iSegment; iVertex<iSegment+2; ++iVertex)
                std::copy_n(
                    triangle.data() + (iVertex % 3) * Dimension,
                    Dimension,
                    segment.data() + (iVertex - iSegment) * Dimension);
            output = output && negativeSide(
                reference,
                std::span<const Scalar,Dimension> (
                    triangle.data() + iSegment * Dimension,
                    Dimension),
                std::span<const Scalar,Dimension> (
                    triangle.data() + ((iSegment + 1) % 3) * Dimension,
                    Dimension));
        }
        return output;
}


} // namespace impl


void MeshData::subdomain(
    std::span<const Scalar> points,
    std::span<DomainData> subdomains) const {
        CIE_CHECK(Dimension * subdomains.size() == points.size(), "")
        CIE_BEGIN_EXCEPTION_TRACING
            for (std::size_t iPoint=0ul; iPoint<subdomains.size(); ++iPoint) {
                const std::span<const Scalar,Dimension> reference(
                    points.data() + iPoint * Dimension,
                    Dimension);
                std::optional<DomainData> maybeDomain;

                for (const auto& [rDomain, rTriangles] : _domainTriangles) {
                    const std::size_t triangleCount = rTriangles.size() / 3 / Dimension;
                    for (std::size_t iTriangle=0ul; iTriangle<triangleCount; ++iTriangle) {
                        const std::span<const Scalar,3*Dimension> triangle(
                            rTriangles.data() + iTriangle * 3 * Dimension,
                            3 * Dimension);
                        if (impl::isInTriangle(reference, triangle)) {
                            maybeDomain = rDomain;
                            break;
                        }
                    } // for iTriangle in triangleCount
                    if (maybeDomain.has_value()) break;
                } // for rDomain, rTriangles in _domainTriangles

                subdomains[iPoint] = maybeDomain.has_value()
                    ? maybeDomain.value()
                    : static_cast<DomainData>(0);
            } // for iPoint in range(subdomains.size())
        CIE_END_EXCEPTION_TRACING
}


std::span<const std::pair<
    MeshData::DomainData,
    Scalar>
> MeshData::domainMap() const noexcept {
    return _domainMap;
}


void MeshData::setCells(RightRef<std::vector<Cell>> rCells) noexcept {
    _cells = std::move(rCells);
}


} // namespace cie::fem
