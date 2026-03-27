// --- Internal Includes ---
#include "poisson2D/BoundaryData.hpp"

// --- STL Includes ---
#include <sstream>


namespace cie::fem {


void io::GraphML::Serializer<BoundaryData>::header(
    Ref<XMLElement> rElement) const {
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        BoundaryData instance;
        this->operator()(defaultElement, instance);
}


void io::GraphML::Serializer<BoundaryData>::operator()(
    Ref<XMLElement> rElement,
    Ref<const BoundaryData> rInstance) const {
        std::ostringstream buffer;
        buffer << rInstance.boundary();
        rElement.addAttribute("bnd", buffer.view());
}


void io::GraphML::Deserializer<BoundaryData>::onElementBegin(
    Ptr<void> pThis,
    std::string_view,
    std::span<AttributePair> attributes) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        const auto it = std::find_if(
            attributes.begin(),
            attributes.end(),
            [] (const auto pair) {
                return pair.first == "bnd";});
        rThis.instance() = BoundaryData(BoundaryID(it->second.data()));
}


void io::GraphML::Deserializer<BoundaryData>::onText(
    Ptr<void>,
    std::string_view) {
        CIE_THROW(
            Exception,
            "Unexpected text block while parsing boundary data.")
}


void io::GraphML::Deserializer<BoundaryData>::onElementEnd(
    Ptr<void> pThis,
    std::string_view elementName) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem
