#ifndef CIE_FEM_GRAPH_HPP
#define CIE_FEM_GRAPH_HPP

// --- External Includes ---
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/stl_extension/inc/NoOpIterator.hpp"
#include "packages/stl_extension/inc/StrongTypeDef.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- STL Includes ---
#include <type_traits> // conditional_t
#include <ranges> // transform_view


namespace cie::fem {



/// @brief A directed graph that automatically manages the connectivities of its
///        @ref Graph::Vertex "vertices" and @ref Graph::Edge "edges".
/// @tparam TVertexData additional data type stored in each vertex.
/// @tparam TEdgeData additional data type stored in each edge.
template <class TVertexData, class TEdgeData>
class Graph
{
private:
    /// @brief Helper class for features common to @ref Vertex and @ref Edge.
    template <class TID>
    class ItemBase
    {
    public:
        ItemBase(TID id) noexcept : _id(id) {}

        TID id() const noexcept {return _id;}

        friend bool operator==(Ref<const ItemBase> rLHS,
                               Ref<const ItemBase> rRHS) noexcept
        {return rLHS.id() == rRHS.id();}

        friend bool operator!=(Ref<const ItemBase> rLHS,
                               Ref<const ItemBase> rRHS) noexcept
        {return rLHS.id() != rRHS.id();}

        struct hash
        {
            std::size_t operator()(Ref<const ItemBase> rItem) const noexcept
            {return std::hash<TID>()(rItem.id());}

            std::size_t operator()(TID id) const noexcept
            {return std::hash<TID>()(id);}
        }; // struct hash

    private:
        ItemBase() = delete;

    private:
        TID _id;
    }; // class ItemBase

public:
    /// @brief Vertex identifier type.
    /// @class VertexID
    /// @see StrongTypeDef
    CIE_STRONG_TYPEDEF(unsigned, VertexID);

    /// @brief Edge identifier type.
    /// @class EdgeID
    /// @see StrongTypeDef
    CIE_STRONG_TYPEDEF(unsigned, EdgeID);

    /// @brief Vertex type of @ref Graph.
    /// @details This vertex implementation stores the IDs of the @ref Edge "edges" originating from-
    ///          or ending at it, as well as additional data associated with it via @p TVertexData.
    ///          Features and storage for the extra data can be disabled without overhead at compile
    ///          time by setting @p TVertexData to @p void.
    class Vertex : public ItemBase<VertexID>
    {
    public:
        /// @copydoc VertexID
        using ID = VertexID;

        /// @brief Construct a new @ref Vertex with an @ref ID and set of @ref Edge "edges".
        /// @param id identifier of the new vertex.
        /// @param rEdges set of edges originating from, or ending at the new vertex.
        Vertex(VertexID id,
               RightRef<tsl::robin_set<EdgeID>> rEdges) noexcept;

        /// @brief Construct a new @ref Vertex with an @ref ID, a set of @ref Edge "edges", and associated data.
        /// @param id identifier of the new vertex.
        /// @param rEdges set of edges originating from, or ending at the new vertex.
        /// @param rData additiona @p VertexData associated with the new vertex.
        /// @note This constructor is disabled if @p TVertexData is @p void.
        Vertex(VertexID id,
               RightRef<tsl::robin_set<EdgeID>> rEdges,
               std::conditional_t<
                    std::is_same_v<std::remove_const_t<TVertexData>,void>,
                    char, // dummy type, cannot be void
                    typename VoidSafe<TVertexData>::RightRef
               > rData) noexcept
        requires (!std::is_same_v<std::remove_const_t<TVertexData>,void>);

        /// @brief Construct a new @ref Vertex with an @ref VertexID, a set of @ref Edge "edges", and associated data.
        /// @param id identifier of the new vertex.
        /// @param rEdges set of edges originating from, or ending at the new vertex.
        /// @param rData additiona @p VertexData associated with the new vertex.
        /// @note This constructor is disabled if @p TVertexData is @p void.
        Vertex(VertexID id,
               RightRef<tsl::robin_set<EdgeID>> rEdges,
               std::conditional_t<
                    std::is_same_v<std::remove_const_t<TVertexData>,void>,
                    int, // dummy type, cannot be void
                    typename VoidSafe<const TVertexData>::Ref
               > rData)
        requires (!std::is_same_v<std::remove_const_t<TVertexData>,void>);

        /// @brief Immutable access to the set of @ref Edge "edge" IDs originating from or ending at this @ref Vertex.
        Ref<const tsl::robin_set<EdgeID>> edges() const noexcept;

        /// @brief Immutable access to the data associated with this @ref Vertex.
        /// @returns an immutable reference to the stored data if @p TVertexData is not @p void,
        ///          otherwise returns @p void.
        typename VoidSafe<const TVertexData>::Ref data() const noexcept;

        /// @brief Mutable access to the data associated with this @ref Vertex.
        /// @returns an immutable reference to the stored data if @p TVertexData is not @p void,
        ///          otherwise returns @p void.
        typename VoidSafe<TVertexData>::Ref data() noexcept;

    private:
        /// @note @ref Graph must be a friend class of @ref Vertex because it needs mutable
        ///            access to vertex connectivities while adding/removing @ref Edge "edges".
        friend class Graph;

        /// @brief Mutable access to the container storing the @ref Edge "edge" IDs associated with this @ref Vertex.
        Ref<tsl::robin_set<EdgeID>> mutableEdges() noexcept;

    private:
        /// @details IDs of @ref Edge "edges" originating from or ending at this @ref Vertex,
        ///          and the additional data associated with it.
        std::conditional_t<
            std::is_same_v<std::remove_const_t<TVertexData>,void>,
            std::tuple<tsl::robin_set<EdgeID>>,
            std::tuple<tsl::robin_set<EdgeID>,TVertexData>
        > _data;
    }; // class Vertex

    /// @brief Edge type of @ref Graph.
    /// @details This edge implementation stores the IDs of its source- and target @ref Vertex "vertices".
    ///          If the edge is directed, it points from the source to the target vertex. Additionally, it
    ///          also stores extra data associated with the edge, which can be disabled without overhead at
    ///          compile time by setting @p TEdgeData to @p void.
    class Edge : public ItemBase<EdgeID>
    {
    public:
        /// @copydoc EdgeID
        using ID = EdgeID;

        /// @brief Construct an @ref Edge from an @ref ID and a source-target @ref Vertex pair.
        /// @param id identifier of the new edge.
        /// @param vertexIDs source and target vertex IDs defining the directed edge.
        Edge(EdgeID id,
             std::pair<VertexID,VertexID> vertexIDs) noexcept;

        /// @brief Construct an @ref Edge from an @ref ID, a source-target @ref Vertex pair, and additional data.
        /// @param id identifier of the new edge.
        /// @param vertexIDs source and target vertex IDs defining the directed edge.
        /// @param rData additional @p TEdgeData associated with the new edge.
        /// @note This constructor is disabled if @p TEdgeData is @p void.
        Edge(EdgeID id,
             std::pair<VertexID,VertexID> vertexIDs,
             std::conditional_t<
                std::is_same_v<std::remove_const_t<TEdgeData>,void>,
                char, // dummy type, cannot be void
                typename VoidSafe<TEdgeData>::RightRef
             > rData) noexcept
        requires (!std::is_same_v<std::remove_const_t<TEdgeData>,void>);

        /// @brief Construct an @ref Edge from an @ref ID, a source-target @ref Vertex pair, and additional data.
        /// @param id identifier of the new edge.
        /// @param vertexIDs source and target vertex IDs defining the directed edge.
        /// @param rData additional @p TEdgeData associated with the new edge.
        /// @note This constructor is disabled if @p TEdgeData is @p void.
        Edge(EdgeID id,
             std::pair<VertexID,VertexID> vertexIDs,
             std::conditional_t<
                std::is_same_v<std::remove_const_t<TEdgeData>,void>,
                int, // dummy type, cannot be void
                typename VoidSafe<const TEdgeData>::Ref
             > rData)
        requires (!std::is_same_v<std::remove_const_t<TEdgeData>,void>);

        /// @brief Immutable access to the source and target @ref Vertex "vertices".
        /// @returns A pair of @ref VertexID "vertex IDs"; first of which belongs to the
        ///          source vertex, while the second one points to the target vertex.
        std::pair<VertexID,VertexID> vertices() const noexcept;

        /// @brief Immutable access to the source @ref Vertex "vertex'" ID.
        VertexID source() const noexcept;

        /// @brief Immutable access to the target @ref Vertex "vertex'" ID.
        VertexID target() const noexcept;

        /// @brief Immutable access to the data associated with this @ref Edge.
        /// @returns an immutable reference to the stored data if @p TEdgeData is not @p void,
        ///          otherwise returns @p void.
        typename VoidSafe<const TEdgeData>::Ref data() const noexcept;

        /// @brief Immutable access to the data associated with this @ref Edge.
        /// @returns an immutable reference to the stored data if @p TEdgeData is not @p void,
        ///          otherwise returns @p void.
        typename VoidSafe<TEdgeData>::Ref data() noexcept;

    private:
        /// @details IDs of the source and target @ref Vertex "vertices", and the additional
        ///          data associated with it.
        std::conditional_t<
            std::is_same_v<std::remove_const_t<TEdgeData>,void>,
            std::tuple<std::pair<VertexID,VertexID>>,
            std::tuple<std::pair<VertexID,VertexID>,TEdgeData>
        > _data;
    }; // class Edge

public:
    /// @brief Construct an empty @ref Graph.
    Graph() noexcept = default;

    /// @brief Move a new @ref Vertex into the @ref Graph.
    /// @param rVertex new vertex to be moved into the graph.
    /// @param overwrite in case of ID clashes, this argument controls whether
    ///                  to replace the existing vertex with the provided one,
    ///                  or do nothing.
    /// @details If a vertex with an identical ID already exists in the graph, two
    ///          things can happen depending on the value of @p overwrite.
    ///          If @p overwrite is @p false (default), then the existing vertex
    ///          remains in the graph and the new one is discarded. If @p overwrite
    ///          is @p true, first the existing vertex is erased (@ref eraseVertex),
    ///          then the new one is inserted.
    /// @see eraseVertex
    /// @returns A mutable reference to the newly inserted, or existing vertex.
    /// @note The new vertex' @ref Edge set must be empty. It will be filled automatically
    ///       by the graph upon edge insertion.
    /// @throws If the edge set of the new vertex is not empty.
    Ref<Vertex> insert(RightRef<Vertex> rVertex,
                       bool overwrite = false);

    /// @brief Copy a new @ref Vertex into the @ref Graph.
    /// @param rVertex new vertex to be copied into the graph.
    /// @param overwrite in case of ID clashes, this argument controls whether
    ///                  to replace the existing vertex with the provided one,
    ///                  or do nothing.
    /// @details If a vertex with an identical ID already exists in the graph, two
    ///          things can happen depending on the value of @p overwrite.
    ///          If @p overwrite is @p false (default), then the existing vertex
    ///          remains in the graph and the new one is discarded. If @p overwrite
    ///          is @p true, first the existing vertex is erased (@ref eraseVertex),
    ///          then the new one is inserted.
    /// @see eraseVertex
    /// @returns A mutable reference to the newly inserted, or existing vertex.
    /// @note The new vertex' @ref Edge set must be empty. It will be filled automatically
    ///       by the graph upon edge insertion.
    /// @throws If the edge set of the new vertex is not empty.
    Ref<Vertex> insert(Ref<const Vertex> rVertex,
                       bool overwrite = false);

    /// @brief Move a new @ref Edge into the @ref Graph.
    /// @param rEdge new edge to be moved into the graph.
    /// @param overwrite in case of ID clashes, this argument controls whether
    ///                  to replace the existing edge with the provided one, or
    ///                  do nothing.
    /// @details If an edge with an identical ID already exists in the graph, two
    ///          things can happen depending on the value of @p overwrite.
    ///          If @p overwrite is @p false (default), then the existing edge
    ///          remains in the graph and the new one is discarded. If @p overwrite
    ///          is @p true, first the existing edge is erased (@ref eraseEdge),
    ///          then the new one is inserted.
    /// @details If the newly inserted edge connects vertices that are not yet part
    ///          of the graph, default constructed ones are inserted.
    /// @see eraseEdge
    Ref<Edge> insert(RightRef<Edge> rEdge,
                     bool overwrite = false);

    /// @brief Copy a new @ref Edge into the @ref Graph.
    /// @param rEdge new edge to be copied into the graph.
    /// @param overwrite in case of ID clashes, this argument controls whether
    ///                  to replace the existing edge with the provided one, or
    ///                  do nothing.
    /// @details If an edge with an identical ID already exists in the graph, two
    ///          things can happen depending on the value of @p overwrite.
    ///          If @p overwrite is @p false (default), then the existing edge
    ///          remains in the graph and the new one is discarded. If @p overwrite
    ///          is @p true, first the existing edge is erased (@ref eraseEdge),
    ///          then the new one is inserted.
    /// @details If the newly inserted edge connects vertices that are not yet part
    ///          of the graph, default constructed ones are inserted.
    /// @see eraseEdge
    Ref<Edge> insert(Ref<const Edge> rEdge,
                     bool overwrite = false);

    /// @brief Erase the @ref Vertex belonging to the provided ID,
    ///        as well as every @ref Edge connected to it.
    /// @param id @ref VertexID "ID" of the vertex to erase.
    /// @details Does nothing if no vertex exists in the @ref Graph
    ///          with the provided @p id.
    /// @see eraseEdge
    /// @returns @p true if a vertex was erased, otherwise @p false.
    bool eraseVertex(VertexID id) noexcept;

    /// @brief Erase the @ref Edge belonging to the provided ID.
    /// @param id @ref EdgeID "ID" of the edge to erase.
    /// @details Also erases the edge ID from the edge set of every
    ///          connected @ref Vertex. Does nothing if no edge exists
    ///          with the provided @p id.
    /// @returns @p true if an edge was erased, otherwise @p false.
    bool eraseEdge(EdgeID id) noexcept;

    /// @brief Find a @ref Vertex in the @ref Graph by its @ref VertexID "ID".
    /// @param id identifier of the vertex to find.
    /// @returns an immutable reference to the vertex with the matching @p id,
    ///          or an empty @ref OptionalRef if no such vertex exists in the
    ///          graph.
    OptionalRef<const Vertex> findVertex(VertexID id) const noexcept;

    /// @brief Find a @ref Vertex in the @ref Graph by its @ref VertexID "ID".
    /// @param id identifier of the vertex to find.
    /// @returns a mutable reference to the vertex with the matching @p id,
    ///          or an empty @ref OptionalRef if no such vertex exists in the
    ///          graph.
    OptionalRef<Vertex> findVertex(VertexID id) noexcept;

    /// @brief Find an @ref Edge in the @ref Graph by its @ref EdgeID "ID".
    /// @param id identifier of the edge to find.
    /// @returns an immutable reference to the edge with the matching @p id,
    ///          or an empty @ref OptionalRef if no such edge exists in the
    ///          graph.
    OptionalRef<const Edge> findEdge(EdgeID id) const noexcept;

    /// @brief Find an @ref Edge in the @ref Graph by its @ref EdgeID "ID".
    /// @param id identifier of the edge to find.
    /// @returns a mutable reference to the edge with the matching @p id,
    ///          or an empty @ref OptionalRef if no such edge exists in the
    ///          graph.
    OptionalRef<Edge> findEdge(EdgeID id) noexcept;

    /// @brief Immutable access to @ref Vertex "vertices" in the @ref Graph.
    /// @returns an immutable view over all vertices.
    auto vertices() const noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_vertices.begin()),
                                  utils::makeNoOpIterator(_vertices.end()),
                                  _vertices.size()),
            [](auto it) -> Ref<const Vertex> {return it.second;}
        );
    }

    /// @brief Mutable access to @ref Vertex "vertices" in the @ref Graph.
    /// @returns a mutable view over all vertices.
    auto vertices() noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_vertices.begin()),
                                  utils::makeNoOpIterator(_vertices.end()),
                                  _vertices.size()),
            [](auto it) -> Ref<Vertex> {return it.value();}
        );
    }

    /// @brief Immutable access to @ref Edge "edges" in the @ref Graph.
    /// @returns an immutable view over all edges.
    auto edges() const noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_edges.begin()),
                                  utils::makeNoOpIterator(_edges.end()),
                                  _edges.size()),
            [](auto it) -> Ref<const Edge> {return it.second;}
        );
    }

    /// @brief Mutable access to @ref Edge "edges" in the @ref Graph.
    /// @returns a mutable view over all vertices.
    auto edges() noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_edges.begin()),
                                  utils::makeNoOpIterator(_edges.end()),
                                  _edges.size()),
            [](auto it) -> Ref<Edge> {return it.value();}
        );
    }

    /// @brief Return true if the graph contains no vertices, false otherwise.
    bool empty() const noexcept;

private:
    tsl::robin_map<
        Size,
        Vertex,
        typename Vertex::hash
    > _vertices;

    tsl::robin_map<
        Size,
        Edge,
        typename Edge::hash
    > _edges;
}; // class Graph


} // namespace cie::fem

#include "packages/graph/impl/Graph_impl.hpp"

#endif
