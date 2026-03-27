#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"
#include "poisson2D/CellData.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/MeshBase.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/io/inc/GraphML.hpp"


namespace cie::fem {


/// @brief Data structure common to the entire @ref Graph "mesh".
class MeshData : public MeshBase<Ansatz> {
public:
    MeshData();

    MeshData(RightRef<Ansatz> rAnsatzSpace);

    CachedQuadraturePointFactory<Dimension,Scalar> makeQuadratureRule() const;

private:
    friend struct io::GraphML::Serializer<MeshData>;

    friend struct io::GraphML::Deserializer<MeshData>;

    /// @brief Set of quadrature points for a default local hypercube.
    /// @details These quadrature points are used while constructing
    ///          cell-specific quadrature rules.
    std::vector<QuadraturePoint<Dimension,Scalar>> _quadraturePointSet;

    DynamicArray<Scalar> _buffer;
}; // class MeshData


/// @brief Serializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Serializer<MeshData> {
    void header(Ref<io::GraphML::XMLElement> rElement) const;

    void operator()(
        [[maybe_unused]] Ref<io::GraphML::XMLElement> rElement,
        [[maybe_unused]] Ref<const MeshData> rInstance) const {}
}; // struct GraphML::Serializer<MeshData>


/// @brief Deserializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Deserializer<MeshData>
    : public io::GraphML::DeserializerBase<MeshData> {
    using io::GraphML::DeserializerBase<MeshData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(
        Ptr<void> pThis,
        std::string_view elementName,
        std::span<io::GraphML::AttributePair>);

    /// @brief This function is called when text block is parsed in the XML document.
    static void onText(
        Ptr<void>,
        std::string_view);

    /// @brief This function is called when an element closing tag is parsed in the XML document.
    static void onElementEnd(
        Ptr<void> pThis,
        std::string_view elementName);

private:
    DynamicArray<Ansatz> _buffer;
}; // struct GraphML::Deserializer<MeshData>


} // namespace cie::fem
