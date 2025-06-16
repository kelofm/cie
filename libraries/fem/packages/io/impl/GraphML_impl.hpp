#ifndef CIE_FEM_GRAPHML_IMPL_HPP
#define CIE_FEM_GRAPHML_IMPL_HPP

// help the language server
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp" // Match::None

// --- STL Includes ---
#include <string> // std::string


namespace cie::fem::io {


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Input::operator()(Ref<Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    SAXHandler sax(this->stream());

    using Deserializer = GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>;
    auto pDeserializer = std::make_unique<Deserializer>();
    sax.push({Deserializer::onElementBegin,
              Deserializer::onText,
              Deserializer::onElementEnd,
              &rGraph});

    sax.parse();

    CIE_END_EXCEPTION_TRACING
}


template <>
struct GraphML::Serializer<void> {};


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Output::operator()(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    XMLElement rootElement = this->root();
    XMLElement graphElement = rootElement.addChild("graph");

    // Write graph attributes.
    graphElement.addAttribute("edgedefault", "directed");

    int dataCount = 0;    // <== counts how many non-void properties need to be written
    std::string graphDataID, vertexDataID, edgeDataID;

    // Write graph data header.
    GraphML::Serializer<TGraphData> graphDataSerializer;

    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        graphDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", graphDataID);
        headerElement.addAttribute("for", "graph");

        graphDataSerializer.header(headerElement);
    }

    // Write vertex data header.
    GraphML::Serializer<TVertexData> vertexDataSerializer;

    if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
        vertexDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", vertexDataID);
        headerElement.addAttribute("for", "node");

        vertexDataSerializer.header(headerElement);
    }

    // Write edge data header.
    GraphML::Serializer<TEdgeData> edgeDataSerializer;

    if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
        edgeDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", edgeDataID);
        headerElement.addAttribute("for", "edge");
        edgeDataSerializer.header(headerElement);
    }

    // Write graph data.
    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        XMLElement dataElement = graphElement.addChild("data");
        dataElement.addAttribute("key", graphDataID);
        graphDataSerializer(dataElement, rGraph.data());
    }

    // Write vertices.
    for (const auto& rItem : rGraph.vertices()) {
        XMLElement element = graphElement.addChild("node");
        element.addAttribute("id", std::to_string(rItem.id()));

        if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            vertexDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", vertexDataID);
        }
    } // for rItem in rGraph.vertices()

    // Write edges.
    for (const auto& rItem : rGraph.edges()) {
        XMLElement element = graphElement.addChild("edge");
        element.addAttribute("id", std::to_string(rItem.id()));
        element.addAttribute("source", std::to_string(rItem.source()));
        element.addAttribute("target", std::to_string(rItem.target()));

        if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            edgeDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", edgeDataID);
        }
    } // for rItem in rGraph.edges()

    // Dump.
    this->write();

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem::io


#endif
