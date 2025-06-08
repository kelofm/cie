// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/io/inc/GraphML.hpp"


namespace cie::fem {


CIE_TEST_CASE("GraphML no data", "[graph]")
{
    CIE_TEST_CASE_INIT("GraphML no data")
    using VertexData = void;
    using EdgeData = void;
    using G = Graph<VertexData,EdgeData>;
    using Edge = G::Edge;

    {
        CIE_TEST_CASE_INIT("output")
        G graph;

        graph.insert(Edge(3, {1, 2}));

        // +---+   3   +---+
        // | 1 |<----->| 2 |
        // +---+       +---+
        io::GraphML::Output io("graphml_no_data_test.graphml");
        CIE_TEST_CHECK_NOTHROW(io(graph));
    }
}


struct GraphMLTestVertexData
{
    std::string data;
}; // struct GraphMLTestVertexData


struct GraphMLTestEdgeData
{
    std::string data;
}; // struct GraphMLTestEdgeData


template <>
struct io::GraphML::Serializer<GraphMLTestVertexData>
{
    void header(Ref<XMLElement> rElement) {
        rElement.addAttribute("id", "GraphMLTestVertexData");
        rElement.addAttribute("for", "node");
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue("default GraphMLTestVertexData");
    }

    void operator()(Ref<XMLElement> rElement, Ref<const GraphMLTestVertexData> rData) {
        rElement.addAttribute("key", "GraphMLTestVertexData");
        rElement.setValue(rData.data);
    }
}; // struct GraphML::Serializer<GraphMLTestVertexData>


template <>
struct io::GraphML::Serializer<GraphMLTestEdgeData>
{
    void header(Ref<XMLElement> rElement) {
        rElement.addAttribute("id", "GraphMLTestEdgeData");
        rElement.addAttribute("for", "edge");
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue("default GraphMLTestEdgeData");
    }

    void operator()(Ref<XMLElement> rElement, Ref<const GraphMLTestEdgeData> rData) {
        rElement.addAttribute("key", "GraphMLTestEdgeData");
        rElement.setValue(rData.data);
    }
}; // struct GraphML::Serializer<GraphMLTestEdgeData>


CIE_TEST_CASE("GraphML custom data", "[graph]")
{
    CIE_TEST_CASE_INIT("GraphML custom data")
    using VertexData = GraphMLTestVertexData;
    using EdgeData = GraphMLTestEdgeData;
    using G = Graph<VertexData,EdgeData>;
    using Edge = G::Edge;

    {
        CIE_TEST_CASE_INIT("output")
        G graph;

        // +---+   1   +---+   2   +---+
        // | 1 |<----->| 2 |<----->| 3 |
        // +---+       +---+       +---+
        graph.insert(Edge(1, {1, 2}));
        graph.insert(Edge(2, {2, 3}));

        const auto rMaybeEdge2 = graph.find(G::EdgeID(2));
        CIE_TEST_CHECK_NOTHROW(rMaybeEdge2.value().data().data = "dummy data of edge 2");

        const auto rMaybeVertex2 = graph.find(G::VertexID(2));
        CIE_TEST_CHECK_NOTHROW(rMaybeVertex2.value().data().data = "dummy data of vertex 2");

        const auto rMaybeVertex3 = graph.find(G::VertexID(3));
        CIE_TEST_CHECK_NOTHROW(rMaybeVertex3.value().data().data = "dummy data of vertex 3");

        io::GraphML::Output io("graphml_custom_data_test.graphml");
        CIE_TEST_CHECK_NOTHROW(io(graph));
    }
}


} // namespace cie::fem
