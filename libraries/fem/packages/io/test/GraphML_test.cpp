// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"


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

    GraphMLTestGraphData(std::string&& rData)
        : data(std::move(rData))
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
struct io::GraphML::Serializer<T>
{
    void header(Ref<XMLElement> rElement) {
        XMLElement defaultData = rElement.addChild("default");
        defaultData.setValue(T().data);
    }

    void operator()(Ref<XMLElement> rElement, Ref<const T> rData) {
        GraphML::XMLElement child = rElement.addChild("value");
        using SubSerializer = GraphML::Serializer<std::string>;
        SubSerializer subSerializer;
        subSerializer(child, rData.data);
    }
}; // GraphML::Serializer<T>


template <class T>
requires (
       std::is_same_v<T,GraphMLTestGraphData>
    || std::is_same_v<T,GraphMLTestVertexData>
    || std::is_same_v<T,GraphMLTestEdgeData>
)
struct io::GraphML::Deserializer<T>
    : public io::GraphML::DeserializerBase<T>
{
    using io::GraphML::DeserializerBase<T>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<AttributePair>) noexcept
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        using SubDeserializer = GraphML::Deserializer<std::string>;
        rThis.sax().push({
            SubDeserializer::make(rThis.instance().data, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    }

    static void onText(Ptr<void>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text block."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName) noexcept
    {
        if (elementName != "value") {
            CIE_THROW(
                Exception,
                "Expecting a closing tag for <value> but got one for <" << elementName << ">."
            )
        }
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }
}; // struct Deserializer<T>


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
        CIE_TEST_REQUIRE_NOTHROW(io(graph));
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


CIE_TEST_CASE("GraphML container", "[graph]")
{
    CIE_TEST_CASE_INIT("GraphML container")

    using G = Graph<
        std::array<std::vector<float>,2>,
        void,
        std::vector<GraphMLTestGraphData>
    >;

    // Test output.
    {
        G graph;

        CIE_TEST_REQUIRE_NOTHROW(
            graph.data() = std::vector<GraphMLTestGraphData> {
                {"dummy data at index 0"},
                {"dummy data at index 1"},
                {"dummy data at index 2"}
            }
        );

        graph.insert(G::Vertex(0)).data() = {
            std::vector<float> {1.2f, 2.4f, 3.6f, 4.8f},
            std::vector<float > {144.0f}};
        graph.insert(G::Vertex(1));
        graph.insert(G::Edge(0, {0, 1}));

        CIE_TEST_CHECK(graph.data()[0].data == "dummy data at index 0");
        CIE_TEST_CHECK(graph.data()[1].data == "dummy data at index 1");
        CIE_TEST_CHECK(graph.data()[2].data == "dummy data at index 2");

        CIE_TEST_REQUIRE_NOTHROW(io::GraphML::Output("graphml_test_vector_data.graphml")(graph));
    }

    // Test input.
    {
        G graph;

        CIE_TEST_CHECK(graph.data().empty());
        CIE_TEST_CHECK(graph.vertices().empty());
        CIE_TEST_CHECK(graph.edges().empty());

        std::ifstream file("graphml_test_vector_data.graphml");
        CIE_TEST_REQUIRE_NOTHROW(io::GraphML::Input(file)(graph));

        CIE_TEST_REQUIRE(graph.data().size() == 3);
        CIE_TEST_CHECK(graph.vertices().size() == 2);
        CIE_TEST_CHECK(graph.edges().size() == 1);

        CIE_TEST_CHECK(graph.data()[0].data == "dummy data at index 0");
        CIE_TEST_CHECK(graph.data()[1].data == "dummy data at index 1");
        CIE_TEST_CHECK(graph.data()[2].data == "dummy data at index 2");

        OptionalRef<G::Vertex> v;

        CIE_TEST_REQUIRE_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_REQUIRE(v.value().data().size() == 2);
        CIE_TEST_REQUIRE(v.value().data()[0].size() == 4);
        CIE_TEST_CHECK(v.value().data()[0][0] == 1.2f);
        CIE_TEST_CHECK(v.value().data()[0][1] == 2.4f);
        CIE_TEST_CHECK(v.value().data()[0][2] == 3.6f);
        CIE_TEST_CHECK(v.value().data()[0][3] == 4.8f);
        CIE_TEST_REQUIRE(v.value().data()[1].size() == 1);
        CIE_TEST_CHECK(v.value().data()[1][0] == 144.0f);

        CIE_TEST_REQUIRE_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().data().size() == 2);
        CIE_TEST_CHECK(v.value().data()[0].empty());
        CIE_TEST_CHECK(v.value().data()[1].empty());

        OptionalRef<G::Edge> e;
        CIE_TEST_REQUIRE_NOTHROW(e = graph.findEdge(0));
        CIE_TEST_CHECK(e.has_value());
    }
}


CIE_TEST_CASE("GraphML heterogeneous", "[graphml]")
{
    CIE_TEST_CASE_INIT("GraphML heterogeneous")

    using Basis = maths::Polynomial<float>;
    using Ansatz = maths::AnsatzSpace<Basis,3>;
    using G = Graph<
        OrientedAxes<3>,
        maths::ScaleTranslateTransform<double,3>,
        Ansatz
    >;

    {
        CIE_TEST_CASE_INIT("output")
        G graph;

        graph.data() = Ansatz({
            Basis({1.0f, 2.5f, 3.8f}),
            Basis({-1.0f, -144.0f})
        });

        std::array<std::array<double,3>,2> transformDefinition {
            std::array<double,3> {1.0, 2.0, 3.0},
            std::array<double,3> {4.0, 5.0, 6.0}
        };

        graph.insert(G::Vertex(0)).data() = OrientedAxes<3>("-z+y-x");
        graph.insert(G::Edge(0, {3, 0})).data() = maths::ScaleTranslateTransform<double,3>(
            transformDefinition.begin(),
            transformDefinition.end()
        );

        CIE_TEST_REQUIRE_NOTHROW(io::GraphML::Output("graphml_test_heterogeneous_data.graphml")(graph));
    }

    {
        CIE_TEST_CASE_INIT("input")
        G graph;

        std::ifstream file("graphml_test_heterogeneous_data.graphml");
        CIE_TEST_REQUIRE_NOTHROW(io::GraphML::Input(file)(graph));

        CIE_TEST_REQUIRE(graph.data().ansatzSet().size() == 2);
        CIE_TEST_REQUIRE(graph.data().ansatzSet()[0].coefficients().size() == 3);
        CIE_TEST_CHECK(graph.data().ansatzSet()[0].coefficients()[0] == 1.0f);
        CIE_TEST_CHECK(graph.data().ansatzSet()[0].coefficients()[1] == 2.5f);
        CIE_TEST_CHECK(graph.data().ansatzSet()[0].coefficients()[2] == 3.8f);
        CIE_TEST_REQUIRE(graph.data().ansatzSet()[1].coefficients().size() == 2);
        CIE_TEST_CHECK(graph.data().ansatzSet()[1].coefficients()[0] == -1.0f);
        CIE_TEST_CHECK(graph.data().ansatzSet()[1].coefficients()[1] == -144.0f);

        OptionalRef<G::Vertex> v;
        CIE_TEST_REQUIRE(graph.vertices().size() == 2);

        CIE_TEST_REQUIRE_NOTHROW(v = graph.find(VertexID(0)));
        CIE_TEST_REQUIRE(v.has_value());
        CIE_TEST_CHECK(v.value().data() == OrientedAxes<3>("-z+y-x"));
        CIE_TEST_CHECK(v.value().data() != OrientedAxes<3>("-z-y-x"));

        CIE_TEST_REQUIRE_NOTHROW(v = graph.find(VertexID(3)));
        CIE_TEST_REQUIRE(v.has_value());
        CIE_TEST_CHECK(v.value().data() == OrientedAxes<3>());

        OptionalRef<G::Edge> e;
        CIE_TEST_REQUIRE(graph.edges().size() == 1);
        CIE_TEST_REQUIRE_NOTHROW(e = graph.find(EdgeID(0)));
        CIE_TEST_REQUIRE(e.has_value());

        std::array<double,3> base       {-1.0, -1.0, -1.0},
                             opposite   { 1.0,  1.0,  1.0};
        std::array<double,3> output;

        CIE_TEST_CHECK_NOTHROW(e.value().data().evaluate(
            base,
            output
        ));
        CIE_TEST_CHECK(output[0] == Catch::Approx(1.0));
        CIE_TEST_CHECK(output[1] == Catch::Approx(2.0));
        CIE_TEST_CHECK(output[2] == Catch::Approx(3.0));

        CIE_TEST_CHECK_NOTHROW(e.value().data().evaluate(
            opposite,
            output
        ));
        CIE_TEST_CHECK(output[0] == Catch::Approx(4.0));
        CIE_TEST_CHECK(output[1] == Catch::Approx(5.0));
        CIE_TEST_CHECK(output[2] == Catch::Approx(6.0));
    }
}


} // namespace cie::fem
