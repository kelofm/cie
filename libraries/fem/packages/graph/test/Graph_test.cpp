// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"


namespace cie::fem {


CIE_TEST_CASE("simple Graph", "[graph]")
{
    CIE_TEST_CASE_INIT("simple Graph")
    using VertexData = void;
    using EdgeData = void;
    using G = Graph<VertexData,EdgeData>;
    using Vertex = G::Vertex;
    using Edge = G::Edge;


    {
        CIE_TEST_CASE_INIT("insert")
        G graph;
        graph.insert(Vertex(1, {}));
        graph.insert(Edge(3, {1, 2}));

        CIE_TEST_CHECK(!graph.findVertex(0).has_value());

        auto r_vertex1 = graph.findVertex(1);
        CIE_TEST_CHECK(r_vertex1.has_value());
        CIE_TEST_CHECK(r_vertex1.value().edges().size() == 1);
        CIE_TEST_CHECK(r_vertex1.value().edges().contains(3));

        auto r_vertex2 = graph.findVertex(2);
        CIE_TEST_CHECK(r_vertex2.has_value());
        CIE_TEST_CHECK(r_vertex2.value().edges().size() == 1);
        CIE_TEST_CHECK(r_vertex2.value().edges().contains(3));

        CIE_TEST_CHECK(!graph.findVertex(3).has_value());

        CIE_TEST_CHECK(!graph.findEdge(0).has_value());
        CIE_TEST_CHECK(!graph.findEdge(1).has_value());
        CIE_TEST_CHECK(!graph.findEdge(2).has_value());

        auto r_edge3 = graph.findEdge(3);
        CIE_TEST_CHECK(r_edge3.has_value());
        CIE_TEST_CHECK(r_edge3.value().source() == 1);
        CIE_TEST_CHECK(r_edge3.value().target() == 2);
    }

    {
        CIE_TEST_CASE_INIT("force insert")
        G graph;
        OptionalRef<Vertex> v;
        OptionalRef<Edge> e;

        graph.insert(Edge(3, {1, 2}));

        // +---+   3   +---+
        // | 1 |<----->| 2 |
        // +---+       +---+
        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(3));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 1);
        CIE_TEST_CHECK(e.value().target() == 2);

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(!v.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(!v.has_value());

        graph.insert(Edge(3, {2, 3})); // <== default behaviour: graph is unchanged

        // +---+   3   +---+
        // | 1 |<----->| 2 |
        // +---+       +---+
        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(3));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 1);
        CIE_TEST_CHECK(e.value().target() == 2);

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(!v.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(!v.has_value());

        graph.insert(Edge(3, {2, 3}), true); // <== force insert: graph is mutated

        // +---+       +---+   3   +---+
        // | 1 |       | 2 |<----->| 3 |
        // +---+       +---+       +---+
        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(3));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 2);
        CIE_TEST_CHECK(e.value().target() == 3);

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(!v.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().empty());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        graph.insert(Vertex(3, {})); // <== default behaviour: graph is unchanged

        // +---+       +---+   3   +---+
        // | 1 |       | 2 |<----->| 3 |
        // +---+       +---+       +---+
        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(3));
        CIE_TEST_CHECK(e.has_value());
        CIE_TEST_CHECK(e.value().source() == 2);
        CIE_TEST_CHECK(e.value().target() == 3);

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(!v.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().empty());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().size() == 1);
        CIE_TEST_CHECK(v.value().edges().contains(3));

        graph.insert(Vertex(3, {}), true); // <== force insert: graph is mutated

        // +---+       +---+
        // | 1 |       | 2 |
        // +---+       +---+
        CIE_TEST_CHECK_NOTHROW(e = graph.findEdge(3));
        CIE_TEST_CHECK(!e.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(0));
        CIE_TEST_CHECK(!v.has_value());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(1));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().empty());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(2));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().empty());

        CIE_TEST_CHECK_NOTHROW(v = graph.findVertex(3));
        CIE_TEST_CHECK(v.has_value());
        CIE_TEST_CHECK(v.value().edges().empty());
    }
}


} // namespace cie::fem
