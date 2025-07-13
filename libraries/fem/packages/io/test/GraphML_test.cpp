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


struct GraphMLTestGraphData
{
    GraphMLTestGraphData()
        : data("default GraphMLTestGraphData")
    {}

    std::string data;
}; // struct GraphMLTestGraphData


struct GraphMLTestVertexData
{
    GraphMLTestVertexData()
        : data("default GraphMLTestVertexData")
    {}


    std::string data;
}; // struct GraphMLTestVertexData


struct GraphMLTestEdgeData
{
    GraphMLTestEdgeData()
        : data("default GraphMLTestEdgeData")
    {}

    std::string data;
}; // struct GraphMLTestEdgeData


template <class T>
requires (
       std::is_same_v<T,GraphMLTestGraphData>
    || std::is_same_v<T,GraphMLTestVertexData>
    || std::is_same_v<T,GraphMLTestEdgeData>
)
struct io::GraphML::Deserializer<T>
{
    using Value = T;

    static void onElementBegin(void*,
                               Ref<SAXHandler>,
                               std::string_view,
                               std::span<AttributePair>) noexcept
    {}

    static void onText(void* pInstance,
                       Ref<SAXHandler>,
                       std::string_view data)
    {
        static_cast<T*>(pInstance)->data = data;
    }

    static void onElementEnd(void*,
                             Ref<SAXHandler>,
                             std::string_view) noexcept
    {}
}; // struct Deserializer<void>


template <>
struct io::GraphML::Serializer<GraphMLTestGraphData>
{
    void header(Ref<XMLElement> rElement) {
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue(GraphMLTestGraphData().data);
    }

    void operator()(Ref<XMLElement> rElement, Ref<const GraphMLTestGraphData> rData) {
        rElement.setValue(rData.data);
    }
}; // GraphML::Serializer<GraphMLTestGraphData>


template <>
struct io::GraphML::Serializer<GraphMLTestVertexData>
{
    void header(Ref<XMLElement> rElement) {
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue(GraphMLTestVertexData().data);
    }

    void operator()(Ref<XMLElement> rElement, Ref<const GraphMLTestVertexData> rData) {
        rElement.setValue(rData.data);
    }
}; // struct GraphML::Serializer<GraphMLTestVertexData>


template <>
struct io::GraphML::Serializer<GraphMLTestEdgeData>
{
    void header(Ref<XMLElement> rElement) {
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue(GraphMLTestEdgeData().data);
    }

    void operator()(Ref<XMLElement> rElement, Ref<const GraphMLTestEdgeData> rData) {
        rElement.setValue(rData.data);
    }
}; // struct GraphML::Serializer<GraphMLTestEdgeData>


CIE_TEST_CASE("GraphML custom data", "[graph]")
{
    CIE_TEST_CASE_INIT("GraphML custom data")
    using GraphData = GraphMLTestGraphData;
    using VertexData = GraphMLTestVertexData;
    using EdgeData = GraphMLTestEdgeData;
    using G = Graph<VertexData,EdgeData,GraphData>;
    using Edge = G::Edge;

    {
        CIE_TEST_CASE_INIT("output")
        G graph;
        graph.data().data = "dummy graph data";

        // +---+   1   +---+   2   +---+
        // | 1 |<----->| 2 |<----->| 3 |
        // +---+       +---+       +---+
        graph.insert(Edge(1, {1, 2}));
        graph.insert(Edge(2, {2, 3}));

        const auto rMaybeEdge2 = graph.find(EdgeID(2));
        CIE_TEST_CHECK_NOTHROW(rMaybeEdge2.value().data().data = "dummy data of edge 2");

        const auto rMaybeVertex2 = graph.find(VertexID(2));
        CIE_TEST_CHECK_NOTHROW(rMaybeVertex2.value().data().data = "dummy data of vertex 2");

        const auto rMaybeVertex3 = graph.find(VertexID(3));
        CIE_TEST_CHECK_NOTHROW(rMaybeVertex3.value().data().data = "dummy data of vertex 3");

        io::GraphML::Output io("graphml_custom_data_test.graphml");
        CIE_TEST_CHECK_NOTHROW(io(graph));
    }

    {
        CIE_TEST_CASE_INIT("input")
        G graph;

        std::ifstream file("graphml_custom_data_test.graphml");
        io::GraphML::Input io(file);

        CIE_TEST_CHECK_NOTHROW(io(graph));
        CIE_TEST_CHECK(graph.data().data == "dummy graph data");

        OptionalRef<G::Vertex> v;
        OptionalRef<G::Edge> e;

        CIE_TEST_CHECK(graph.vertices().size() == 3);

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().data().data == "default GraphMLTestVertexData");
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().find(EdgeID(1)) != v.value().edges().end());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().data().data == "dummy data of vertex 2");
        CIE_TEST_CHECK(v.value().edges().size() == 2);
        CIE_TEST_CHECK(v.value().edges().find(EdgeID(1)) != v.value().edges().end());
        CIE_TEST_CHECK(v.value().edges().find(EdgeID(2)) != v.value().edges().end());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().data().data == "dummy data of vertex 3");
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().find(EdgeID(2)) != v.value().edges().end());

        CIE_TEST_CHECK(graph.edges().size() == 2);

        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(EdgeID(1)));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 1);
        CIE_TEST_CHECK(e.value().target() == 2);
        CIE_TEST_CHECK(e.value().data().data == "default GraphMLTestEdgeData");

        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(EdgeID(2)));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 2);
        CIE_TEST_CHECK(e.value().target() == 3);
        CIE_TEST_CHECK(e.value().data().data == "dummy data of edge 2");
    }
}


} // namespace cie::fem
