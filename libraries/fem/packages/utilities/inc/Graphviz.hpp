#ifndef CIE_FEM_IO_GRAPHVIZ_HPP
#define CIE_FEM_IO_GRAPHVIZ_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"

// --- Utiltiy Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <iostream> // istream, ostream


namespace cie::io {


/// @brief IO class for reading and writing matrices in @a DOT format.
/// @ingroup cieutils
struct Graphviz
{
    struct Settings
    {
        bool directed = true;
    }; // struct Settings

    class Input
    {}; // class Input

    class Output
    {
    public:
        /// @brief Construct a DOT output object writing to @p stdout.
        Output();

        /// @brief Constrtuct a DOT output object writing to the provided stream.
        Output(Ref<std::ostream> rStream, Settings settings = {.directed = true});

        template <class TVertexData, class TEdgeData>
        Ref<Output> operator()(Ref<const fem::Graph<TVertexData,TEdgeData>> rGraph);

        template <class TVertexData, class TEdgeData, class TVertexFunctor, class TEdgeFunctor>
        Ref<Output> operator()(Ref<const fem::Graph<TVertexData,TEdgeData>> rGraph,
                               TVertexFunctor&& rVertexFunctor,
                               TEdgeFunctor&& rEdgeFunctor);

    private:
        Ptr<std::ostream> _pStream;

        Settings _settings;
    }; // class Output
}; // struct Graphviz


} // namespace cie::io

#include "packages/utilities/impl/Graphviz_impl.hpp"

#endif
