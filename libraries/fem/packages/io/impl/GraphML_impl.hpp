#ifndef CIE_FEM_GRAPHML_IMPL_HPP
#define CIE_FEM_GRAPHML_IMPL_HPP

// help the language server
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp" // Match::None


namespace cie::fem::io {


template <>
struct GraphML::Serializer<void> {};


template <class TVertexData, class TEdgeData>
void GraphML::Output::operator()(Ref<const Graph<TVertexData,TEdgeData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    XMLElement rootElement = this->root();
    XMLElement graphElement = rootElement.addChild("graph");

    // Write graph attributes.
    graphElement.addAttribute("edgedefault", "directed");

    // Write vertex data header.
    GraphML::Serializer<TVertexData> vertexDataSerializer;

    if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
        XMLElement vertexHeader = graphElement.addChild("key");
        vertexHeader.addAttribute("id", "0");
        vertexHeader.addAttribute("for", "node");
        vertexDataSerializer.header(vertexHeader);
    }

    // Write edge data header.
    GraphML::Serializer<TEdgeData> edgeDataSerializer;

    if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
        XMLElement edgeHeader = graphElement.addChild("key");
        edgeHeader.addAttribute("id", "1");
        edgeHeader.addAttribute("for", "edge");
        edgeDataSerializer.header(edgeHeader);
    }

    // Write vertices.
    for (const auto& rVertex : rGraph.vertices()) {
        XMLElement vertexElement = graphElement.addChild("node");
        vertexElement.addAttribute("id", std::to_string(rVertex.id()));

        if constexpr (!std::is_same_v<TVertexData,void>) {
            XMLElement vertexDataElement = vertexElement.addChild("data");
            vertexDataSerializer(vertexDataElement, rVertex.data());
            vertexDataElement.addAttribute("key", "0");
        }
    } // for rVertex in rGraph.vertices()

    // Write edges.
    for (const auto& rEdge : rGraph.edges()) {
        XMLElement edgeElement = graphElement.addChild("edge");
        edgeElement.addAttribute("id", std::to_string(rEdge.id()));
        edgeElement.addAttribute("source", std::to_string(rEdge.source()));
        edgeElement.addAttribute("target", std::to_string(rEdge.target()));

        if constexpr (!std::is_same_v<TEdgeData,void>) {
            XMLElement edgeDataElement = edgeElement.addChild("data");
            edgeDataSerializer(edgeDataElement, rEdge.data());
            edgeDataElement.addAttribute("key", "1");
        }
    } // for rEdge in rGraph.edges()

    // Dump.
    this->write();

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem::io


#endif
