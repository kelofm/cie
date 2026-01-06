#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"
#include "poisson2D/CellData.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"


namespace cie::fem {


/// @brief Data structure common to the entire @ref Graph "mesh".
class MeshData {
public:
    MeshData() noexcept = default;

    MeshData(RightRef<DynamicArray<Ansatz>> rAnsatzSpaces,
             std::span<unsigned> integrationOrders)
        : _ansatzSpaces(std::move(rAnsatzSpaces)),
          _ansatzDerivatives(),
          _quadraturePointSets()
    {
        // Cache the derivatives of ansatz spaces.
        _ansatzDerivatives.reserve(_ansatzSpaces.size());
        std::transform(
            _ansatzSpaces.begin(),
            _ansatzSpaces.end(),
            std::back_inserter(_ansatzDerivatives),
            [](Ref<const Ansatz> rAnsatz){return rAnsatz.makeDerivative();});

        // Cache quadrature points.
        _quadraturePointSets.resize(integrationOrders.size());
        CellData dummyElement;
        for (std::size_t iSet=0u; iSet<integrationOrders.size(); ++iSet) {
            const std::size_t integrationOrder = integrationOrders[iSet];

            // Generate 1D quadrature points.
            DynamicArray<QuadraturePoint<1,Scalar>> basePoints;
            OuterProductQuadraturePointFactory<Dimension,Scalar> generator;

            {
                GaussLegendreQuadrature<Scalar> quadrature(integrationOrder);
                basePoints.reserve(quadrature.numberOfNodes());
                for (std::size_t iNode=0ul; iNode<quadrature.numberOfNodes(); ++iNode) {
                    basePoints.emplace_back(
                        quadrature.nodes()[iNode],
                        quadrature.weights()[iNode]);
                }
                generator = OuterProductQuadraturePointFactory<Dimension,Scalar>(basePoints);
            }

            // Generate nD quadrature points.
            constexpr std::size_t initialSize = 0x10;
            constexpr double growFactor = 2.0;
            std::size_t pointCount = 0ul;
            auto& rSet = _quadraturePointSets[iSet];
            rSet.resize(initialSize);

            while (true) {
                // Extend the container if necessary.
                if (!(pointCount < rSet.size())) {
                    rSet.resize(growFactor * rSet.size());
                }

                std::span<QuadraturePoint<Dimension,Scalar>> targetSpan(
                    rSet.data() + pointCount,
                    rSet.data() + rSet.size());

                // Request a new batch of quadrature points.
                const auto newPointCount = generator.generate(dummyElement, targetSpan);
                pointCount += newPointCount;

                if (!newPointCount) break;
            }

            rSet.resize(pointCount);
        } // for iSet in range(integrationOrders.size())
    }

    std::span<const Ansatz> ansatzSpaces() const noexcept {
        return _ansatzSpaces;
    }

    std::span<const Ansatz::Derivative> ansatzDerivatives() const noexcept {
        return _ansatzDerivatives;
    }

    CachedQuadraturePointFactory<Dimension,Scalar> makeQuadratureRule(std::size_t iAnsatz) const {
        return CachedQuadraturePointFactory<Dimension,Scalar>(
            std::span<const QuadraturePoint<Dimension,Scalar>>(
                _quadraturePointSets[iAnsatz].data(),
                _quadraturePointSets[iAnsatz].size())
        );
    }

private:
    friend struct io::GraphML::Serializer<MeshData>;

    friend struct io::GraphML::Deserializer<MeshData>;

    /// @brief Collection of all ansatz spaces the contained cells can refer to.
    DynamicArray<Ansatz> _ansatzSpaces;

    /// @brief Collection of all ansatz spaces' derivatives the contained cells can refer to.
    DynamicArray<Ansatz::Derivative> _ansatzDerivatives;

    DynamicArray<DynamicArray<QuadraturePoint<Dimension,Scalar>>> _quadraturePointSets;
}; // class MeshData


/// @brief Serializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Serializer<MeshData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const {
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

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const MeshData> rInstance) const {
        io::GraphML::Serializer<DynamicArray<Ansatz>> subSerializer;
        subSerializer(rElement, rInstance._ansatzSpaces);
    }
}; // struct GraphML::Serializer<MeshData>

/// @brief Deserializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Deserializer<MeshData>
    : public io::GraphML::DeserializerBase<MeshData>
{
    using io::GraphML::DeserializerBase<MeshData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair>) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        // Defer parsing to the array of ansatz spaces.
        using SubDeserializer = io::GraphML::Deserializer<DynamicArray<Ansatz>>;
        rThis.sax().push({
            SubDeserializer::make(rThis._buffer, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    }

    /// @brief This function is called when text block is parsed in the XML document.
    static void onText(void*,
                       std::string_view) {
        // No text data is expected for this class.
        CIE_THROW(Exception, "Unexpected text block while parsing mesh data.")
    }

    /// @brief This function is called when an element closing tag is parsed in the XML document.
    static void onElementEnd(void* pThis,
                             std::string_view elementName) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        // Move the parsed ansatz spaces from the buffer to the mesh data instance.
        rThis.instance()._ansatzSpaces = std::move(rThis._buffer);

        // Build derivatives from the parsed ansatz spaces.
        std::transform(rThis.instance()._ansatzSpaces.begin(),
                       rThis.instance()._ansatzSpaces.end(),
                       std::back_inserter(rThis.instance()._ansatzDerivatives),
                       [] (Ref<const Ansatz> rAnsatz) -> Ansatz::Derivative {
                            return rAnsatz.makeDerivative();
                       });

        // The parser's job is done => destroy it.
        rThis.template release<Deserializer>(&rThis, elementName);
    }

private:
    DynamicArray<Ansatz> _buffer;
}; // struct GraphML::Deserializer<MeshData>


} // namespace cie::fem
