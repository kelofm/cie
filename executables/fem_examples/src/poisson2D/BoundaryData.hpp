#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"


namespace cie::fem {


/// @brief Data structure unique to @ref Graph::Edge "boundaries" between cells.
struct BoundaryData {
    static constexpr inline unsigned Dimension = 2;

    constexpr BoundaryData() noexcept = default;

    constexpr BoundaryData(BoundaryID boundary) noexcept
        : _boundary(boundary)
    {}

    /// @brief Boundary identifier of the shared boundary between the adjacent cells.
    constexpr BoundaryID boundary() const noexcept {
        return _boundary;
    }

private:
    BoundaryID _boundary;
}; // struct BoundaryData


template <>
struct io::GraphML::Serializer<BoundaryData> {
    void header(Ref<io::GraphML::XMLElement> rElement) const;

    void operator()(
        Ref<io::GraphML::XMLElement> rElement,
        Ref<const BoundaryData> rInstance) const;
}; // struct GraphML::Serializer<BoundaryData>


template <>
struct io::GraphML::Deserializer<BoundaryData>
    : public io::GraphML::DeserializerBase<BoundaryData> {
    using io::GraphML::DeserializerBase<BoundaryData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(
        void* pThis,
        std::string_view,
        std::span<io::GraphML::AttributePair> attributes);

    static void onText(Ptr<void>, std::string_view);

    static void onElementEnd(
        Ptr<void> pThis,
        std::string_view elementName);
}; // struct GraphML::Deserializer<BoundaryData>


} // namespace cie::fem
