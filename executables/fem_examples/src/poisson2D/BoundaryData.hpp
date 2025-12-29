#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"


namespace cie::fem {


/// @brief Data structure unique to @ref Graph::Edge "boundaries" between cells.
struct BoundaryData
{
    /// @brief Boundary identifier of the shared boundary between the adjacent cells.
    BoundaryID boundary;
}; // struct BoundaryData


template <>
struct io::GraphML::Serializer<BoundaryData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const {
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        BoundaryData instance;
        this->operator()(defaultElement, instance);
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const BoundaryData> rInstance) const {
        std::ostringstream buffer;
        buffer << rInstance.boundary;
        rElement.addAttribute("bnd", buffer.view());
    }
}; // struct GraphML::Serializer<BoundaryData>


template <>
struct io::GraphML::Deserializer<BoundaryData>
    : public io::GraphML::DeserializerBase<BoundaryData>
{
    using io::GraphML::DeserializerBase<BoundaryData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view,
                               std::span<io::GraphML::AttributePair> attributes) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                            return pair.first == "bnd";
                                         });
            rThis.instance().boundary = BoundaryID(it->second.data());
        }
    }

    static void onText(Ptr<void>, std::string_view) {
        CIE_THROW(
            Exception,
            "Unexpected text block while parsing boundary data."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }

}; // struct GraphML::Deserializer<BoundaryData>


} // namespace cie::fem
