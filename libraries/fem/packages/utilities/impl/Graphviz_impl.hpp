#ifndef CIE_FEM_IO_GRAPHVIZ_IMPL_HPP
#define CIE_FEM_IO_GRAPHVIZ_IMPL_HPP

// Help the language server.
#include "packages/utilities/inc/Graphviz.hpp"

// --- STL Includes ---
#include <string>


namespace cie::io {


template <class TVertexData, class TEdgeData>
Ref<Graphviz::Output> Graphviz::Output::operator()(Ref<const fem::Graph<TVertexData,TEdgeData>> rGraph)
{
    return this->operator()(rGraph,
                            [](const auto&){return "";},
                            [](const auto&){return "";});
}


template <class TVertexData, class TEdgeData, class TVertexFunctor, class TEdgeFunctor>
Ref<Graphviz::Output> Graphviz::Output::operator()(Ref<const fem::Graph<TVertexData,TEdgeData>> rGraph,
                                                   TVertexFunctor&& rVertexFunctor,
                                                   TEdgeFunctor&& rEdgeFunctor)
{
    Ref<std::ostream> rStream = *_pStream;
    const std::string indent("    ");
    const std::string edge = _settings.directed ? "->" : "--";

    // Write header
    rStream << (_settings.directed ? "digraph" : "graph") << " graphname {\n";

    // Write vertices
    for (const auto& rVertex : rGraph.vertices()) {
        std::string label;

        CIE_BEGIN_EXCEPTION_TRACING
        label = rVertexFunctor(rVertex);
        CIE_END_EXCEPTION_TRACING

        rStream << indent
                << rVertex.id();

        if (!label.empty()) {
            rStream << " [label=\""
                    << label
                    << "\"]";
        }

        rStream << ";\n";
    } // for rVertex in rGraph.vertices

    // Write edges
    for (const auto& rEdge : rGraph.edges()) {
        std::string label;

        CIE_BEGIN_EXCEPTION_TRACING
        label = rEdgeFunctor(rEdge);
        CIE_END_EXCEPTION_TRACING

        rStream << indent
                << rEdge.source() << ' ' << edge << ' ' << rEdge.target();

        if (!label.empty()) {
            rStream << " [label=\""
                    << label
                    << "\"]";
        }

        rStream << ";\n";
    } // for rEdge in rGraph.edges

    // Write footer
    rStream << "}\n";

    return *this;
}


} // namespace cie::io


#endif
