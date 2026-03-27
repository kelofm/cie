// --- Internal Includes ---
#include "poisson2D/MeshData.hpp"

// --- FEM Includes ---
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"


namespace cie::fem {


MeshData::MeshData()
    : MeshBase<Ansatz>(),
      _quadraturePointSet(),
      _buffer()
{}


MeshData::MeshData(RightRef<Ansatz> rAnsatzSpace)
    : MeshBase<Ansatz>(std::span<const Ansatz>(&rAnsatzSpace, 1)),
      _quadraturePointSet(),
      _buffer() {
        // Cache quadrature points.
        // Generate 1D quadrature points.
        DynamicArray<QuadraturePoint<1,Scalar>> basePoints;
        OuterProductQuadraturePointFactory<Dimension,Scalar> generator;

        CIE_BEGIN_EXCEPTION_TRACING
            GaussLegendreQuadrature<Scalar> quadrature(integrationOrder);
            basePoints.reserve(quadrature.numberOfNodes());
            for (std::size_t iNode=0ul; iNode<quadrature.numberOfNodes(); ++iNode) {
                basePoints.emplace_back(
                    quadrature.nodes()[iNode],
                    quadrature.weights()[iNode]);
            }
            generator = OuterProductQuadraturePointFactory<Dimension,Scalar>(basePoints);
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

            std::span<QuadraturePoint<Dimension,Scalar>> targetSpan(
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


CachedQuadraturePointFactory<Dimension,Scalar> MeshData::makeQuadratureRule() const {
    return CachedQuadraturePointFactory<Dimension,Scalar>(
        std::span<const QuadraturePoint<Dimension,Scalar>>(
            _quadraturePointSet.data(),
            _quadraturePointSet.size())
    );
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
