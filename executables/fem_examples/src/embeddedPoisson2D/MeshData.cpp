// --- Internal Includes ---
#include "embeddedPoisson2D/MeshData.hpp"

// --- FEM Includes ---
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include <limits>


namespace cie::fem {


MeshData::MeshData()
    : MeshBase<Ansatz>(),
      _quadraturePointSet()
{}


MeshData::MeshData(
    RightRef<Ansatz> rAnsatzSpace,
    RightRef<std::vector<std::pair<DomainData,std::vector<Scalar>>>> domainTriangles,
    std::span<const std::pair<DomainData,Scalar>> domainMap)
        :   MeshBase<Ansatz>(std::span<const Ansatz>(&rAnsatzSpace, 1)),
            _quadraturePointSet(),
            _domainTriangles(std::move(domainTriangles)),
            _domainMap(domainMap.begin(), domainMap.end()) {
                // Cache quadrature points.
                // Generate 1D quadrature points.
                DynamicArray<QuadraturePoint<1,Scalar,Scalar>> basePoints;
                OuterProductQuadraturePointFactory<Dimension,Scalar,Scalar> generator;

                CIE_BEGIN_EXCEPTION_TRACING
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
    CellData,
    Scalar
> MeshData::makeQuadratureRule(Ref<const CellData> rCell) const {
    return KDTreeQuadraturePointFactory<
        Dimension,
        Scalar,
        CellData,
        Scalar>(
            *this,
            rCell,
            std::span<const QuadraturePoint<Dimension,Scalar,Scalar>>(
                _quadraturePointSet.data(),
                _quadraturePointSet.size()),
        {0ul, 6ul}
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


void io::GraphML::Serializer<MeshData>::header(Ref<io::GraphML::XMLElement> rElement) const {
    // Add default value to the header.
    io::GraphML::XMLElement defaultElement = rElement.addChild("default");
    MeshData instance;
    this->operator()(defaultElement, instance);

    // Add description to the header.
    io::GraphML::XMLElement descriptionElement = rElement.addChild("desc");
    std::stringstream description;
    description << "Data structure shared by all cells and boundaries of the mesh. "
                << "In this case, this means the ansatz spaces of the cells as "
                << "well as their derivatives. Each cell stores an index referring "
                << "to their own ansatz spaces.",
    descriptionElement.setValue(description.view());
}


void io::GraphML::Deserializer<MeshData>::onElementBegin(
    Ptr<void> pThis,
    std::string_view elementName,
    std::span<io::GraphML::AttributePair>) {
        // TODO
        (void)(pThis);
        CIE_THROW(NotImplementedException, elementName)
//        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
//
//        // Defer parsing to the array of ansatz spaces.
//        using SubDeserializer = io::GraphML::Deserializer<DynamicArray<Ansatz>>;
//        rThis.sax().push({
//            SubDeserializer::make(rThis._buffer, rThis.sax(), elementName),
//            SubDeserializer::onElementBegin,
//            SubDeserializer::onText,
//            SubDeserializer::onElementEnd
//        });
}


void io::GraphML::Deserializer<MeshData>::onText(
    Ptr<void>,
    std::string_view) {
        // No text data is expected for this class.
        CIE_THROW(Exception, "Unexpected text block while parsing mesh data.")
}


void io::GraphML::Deserializer<MeshData>::onElementEnd(
    Ptr<void> pThis,
    std::string_view elementName) {
        // TODO
        (void)(pThis);
        CIE_THROW(NotImplementedException, elementName)
//        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
//
//        // Move the parsed ansatz spaces from the buffer to the mesh data instance.
//        rThis.instance()._ansatzSpaces = std::move(rThis._buffer);
//
//        // Build derivatives from the parsed ansatz spaces.
//        std::transform(rThis.instance()._ansatzSpaces.begin(),
//                       rThis.instance()._ansatzSpaces.end(),
//                       std::back_inserter(rThis.instance()._ansatzDerivatives),
//                       [] (Ref<const Ansatz> rAnsatz) -> AnsatzDerivative {
//                            return rAnsatz.makeDerivative();
//                       });
//
//        // The parser's job is done => destroy it.
//        rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem
