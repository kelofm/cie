#ifndef CIE_FEM_GRAPH_HPP
#define CIE_FEM_GRAPH_HPP

// --- External Includes ---
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/stl_extension/inc/NoOpIterator.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- STL Includes ---
#include <type_traits> // conditional_t
#include <ranges> // transform_view


namespace cie::fem {



template <class TVertexData, class TEdgeData>
class Graph
{
private:
    class ItemBase
    {
    public:
        ItemBase(Size id) noexcept : _id(id) {}

        Size id() const noexcept {return _id;}

        friend bool operator==(Ref<const ItemBase> r_lhs,
                               Ref<const ItemBase> r_rhs) noexcept
        {return r_lhs.id() == r_rhs.id();}

        friend bool operator!=(Ref<const ItemBase> r_lhs,
                               Ref<const ItemBase> r_rhs) noexcept
        {return r_lhs.id() != r_rhs.id();}

        struct hash
        {
            std::size_t operator()(Ref<const ItemBase> r_item) const noexcept
            {return r_item.id();}
        }; // struct hash

    private:
        ItemBase() = delete;

    private:
        Size _id;
    }; // class ItemBase

public:
    struct Vertex : public ItemBase
    {
        Vertex(Size id,
               RightRef<tsl::robin_set<Size>> r_edges) noexcept;

        Vertex(Size id,
               RightRef<tsl::robin_set<Size>> r_edges,
               std::conditional_t<
                    std::is_same_v<std::remove_const_t<TVertexData>,void>,
                    char, // dummy type, cannot be void
                    typename VoidSafe<TVertexData>::RightRef
               > r_data) noexcept
        requires (!std::is_same_v<std::remove_const_t<TVertexData>,void>);

        Vertex(Size id,
               RightRef<tsl::robin_set<Size>> r_edges,
               std::conditional_t<
                    std::is_same_v<std::remove_const_t<TVertexData>,void>,
                    int, // dummy type, cannot be void
                    typename VoidSafe<const TVertexData>::Ref
               > r_data)
        requires (!std::is_same_v<std::remove_const_t<TVertexData>,void>);

        Ref<const tsl::robin_set<Size>> edges() const noexcept;

        typename VoidSafe<const TVertexData>::Ref data() const noexcept;

        typename VoidSafe<TVertexData>::Ref data() noexcept;

    private:
        friend class Graph;

        Ref<tsl::robin_set<Size>> mutableEdges() noexcept;

    private:
        std::conditional_t<
            std::is_same_v<std::remove_const_t<TVertexData>,void>,
            std::tuple<tsl::robin_set<Size>>,
            std::tuple<tsl::robin_set<Size>,TVertexData>
        > _data;
    }; // class Vertex

    class Edge : public ItemBase
    {
    public:
        Edge(Size id,
             std::pair<Size,Size> vertexIDs) noexcept;

        Edge(Size id,
             std::pair<Size,Size> vertexIDs,
             std::conditional_t<
                std::is_same_v<std::remove_const_t<TEdgeData>,void>,
                char, // dummy type, cannot be void
                typename VoidSafe<TEdgeData>::RightRef
             > r_data) noexcept
        requires (!std::is_same_v<std::remove_const_t<TEdgeData>,void>);

        Edge(Size id,
             std::pair<Size,Size> vertexIDs,
             std::conditional_t<
                std::is_same_v<std::remove_const_t<TEdgeData>,void>,
                int, // dummy type, cannot be void
                typename VoidSafe<const TEdgeData>::Ref
             > r_data)
        requires (!std::is_same_v<std::remove_const_t<TEdgeData>,void>);

        std::pair<Size,Size> vertices() const noexcept;

        Size source() const noexcept;

        Size target() const noexcept;

        typename VoidSafe<const TEdgeData>::Ref data() const noexcept;

        typename VoidSafe<TEdgeData>::Ref data() noexcept;

    private:
        std::conditional_t<
            std::is_same_v<std::remove_const_t<TEdgeData>,void>,
            std::tuple<std::pair<Size,Size>>,
            std::tuple<std::pair<Size,Size>,TEdgeData>
        > _data;
    }; // class Edge

public:
    Graph() noexcept = default;

    Ref<Vertex> insert(Ref<const Vertex> r_data,
                       bool overwrite = false);

    Ref<Edge> insert(Ref<const Edge> r_data,
                     bool overwrite = false);

    bool eraseVertex(Size id) noexcept;

    bool eraseEdge(Size id) noexcept;

    OptionalRef<const Vertex> findVertex(Size id) const noexcept;

    OptionalRef<Vertex> findVertex(Size id) noexcept;

    OptionalRef<const Edge> findEdge(Size id) const noexcept;

    OptionalRef<Edge> findEdge(Size id) noexcept;

    auto vertices() const noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_vertices.begin()),
                                  utils::makeNoOpIterator(_vertices.end())),
            [](auto it) -> Ref<const Vertex> {return it.second;}
        );
    }

    auto vertices() noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_vertices.begin()),
                                  utils::makeNoOpIterator(_vertices.end())),
            [](auto it) -> Ref<Vertex> {return it.value();}
        );
    }

    auto edges() const noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_edges.begin()),
                                  utils::makeNoOpIterator(_edges.end())),
            [](auto it) -> Ref<const Edge> {return it.second;}
        );
    }

    auto edges() noexcept
    {
        return std::ranges::transform_view(
            std::ranges::subrange(utils::makeNoOpIterator(_edges.begin()),
                                  utils::makeNoOpIterator(_edges.end())),
            [](auto it) -> Ref<Edge> {return it.value();}
        );
    }

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
