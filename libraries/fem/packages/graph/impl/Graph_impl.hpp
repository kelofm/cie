#ifndef CIE_FEM_GRAPH_IMPL_HPP
#define CIE_FEM_GRAPH_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/types/inc/modifiers.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <type_traits>


namespace cie::fem {


namespace impl {


template <class TContainer>
OptionalRef<typename CopyConstQualifier<
    TContainer,
    typename std::remove_reference_t<TContainer>::mapped_type
>::Type> findGraphItem(Size id, TContainer&& r_container) noexcept
{
    const auto it = r_container.find(id);
    if (it == r_container.end()) {
        return {};
    } else {
        return it.value();
    }
}


} // namespace impl


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(Size id,
                               RightRef<tsl::robin_set<Size>> r_edges) noexcept
    : ItemBase(id),
      _data(std::move(r_edges))
{}


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(Size id,
                               RightRef<tsl::robin_set<Size>> r_edges,
                               std::conditional_t<
                                    std::is_same_v<std::remove_const_t<TVD>,void>,
                                    char, // dummy type, cannot be void
                                    typename VoidSafe<TVD>::RightRef
                               > r_data) noexcept
requires (!std::is_same_v<std::remove_const_t<TVD>,void>)
    : ItemBase(id),
      _data(std::move(r_edges), std::move(r_data))
{}


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(Size id,
                               RightRef<tsl::robin_set<Size>> r_edges,
                               std::conditional_t<
                                    std::is_same_v<std::remove_const_t<TVD>,void>,
                                    int, // dummy type, cannot be void
                                    typename VoidSafe<const TVD>::Ref
                               > r_data)
requires (!std::is_same_v<std::remove_const_t<TVD>,void>)
    : ItemBase(id),
      _data(std::move(r_edges), r_data)
{}


template <class TVD, class TED>
Ref<const tsl::robin_set<Size>> Graph<TVD,TED>::Vertex::edges() const noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
typename VoidSafe<const TVD>::Ref Graph<TVD,TED>::Vertex::data() const noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TVD>,void>) {
        std::get<1>(_data);
    }
}


template <class TVD, class TED>
typename VoidSafe<TVD>::Ref Graph<TVD,TED>::Vertex::data() noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TVD>,void>) {
        std::get<1>(_data);
    }
}


template <class TVD, class TED>
Ref<tsl::robin_set<Size>> Graph<TVD,TED>::Vertex::mutableEdges() noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(Size id,
                           std::pair<Size,Size> vertices) noexcept
    : ItemBase(id),
      _data(vertices)
{
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(Size id,
                           std::pair<Size,Size> vertices,
                           std::conditional_t<
                                std::is_same_v<std::remove_const_t<TED>,void>,
                                char, // dummy type, cannot be void
                                typename VoidSafe<TED>::RightRef
                           > r_data) noexcept
requires (!std::is_same_v<std::remove_const_t<TED>,void>)
    : ItemBase(id),
      _data(vertices, std::move(r_data))
{
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(Size id,
                           std::pair<Size,Size> vertices,
                           std::conditional_t<
                                std::is_same_v<std::remove_const_t<TED>,void>,
                                int, // dummy type, cannot be void
                                typename VoidSafe<const TED>::Ref
                           > r_data)
requires (!std::is_same_v<std::remove_const_t<TED>,void>)
    : ItemBase(id),
      _data(vertices, r_data)
{
}


template <class TVD, class TED>
std::pair<Size,Size> Graph<TVD,TED>::Edge::vertices() const noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
Size Graph<TVD,TED>::Edge::source() const noexcept
{
    return std::get<0>(_data).first;
}


template <class TVD, class TED>
Size Graph<TVD,TED>::Edge::target() const noexcept
{
    return std::get<0>(_data).second;
}


template <class TVD, class TED>
typename VoidSafe<const TED>::Ref Graph<TVD,TED>::Edge::data() const noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TED>,void>) {
        return std::get<1>(_data);
    }
}


template <class TVD, class TED>
typename VoidSafe<TED>::Ref Graph<TVD,TED>::Edge::data() noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TED>,void>) {
        return std::get<1>(_data);
    }
}


template <class TVD, class TED>
Ref<typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::insert(Ref<const Vertex> r_vertex,
                                                            bool overwrite)
{
    CIE_CHECK(r_vertex.edges().empty(), "Attempt to insert vertex " << r_vertex.id() << " that already has edges")

    const auto id = r_vertex.id();
    std::pair<typename decltype(_vertices)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no vertex exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseVertex(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _vertices.emplace(id, r_vertex);
    CIE_END_EXCEPTION_TRACING

    // Check whether every edge connected to
    // the newly inserted vertex exists.
    if (emplaceResult.second) {
        for (Size edgeID : emplaceResult.first.value().edges()) {
            CIE_CHECK(
                this->findEdge(edgeID).has_value(),
                "Vertex " << id << " is connected to edge " << edgeID << " that is not part of the graph\n"
            );
        }
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
Ref<typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::insert(Ref<const Edge> r_edge,
                                                          bool overwrite)
{
    const auto id = r_edge.id();
    std::pair<typename decltype(_edges)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no edge exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseEdge(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _edges.emplace(id, r_edge);
    CIE_END_EXCEPTION_TRACING

    // Ensure that both endpoints of the edge exist in the graph,
    // and both end points reference the edge.
    if (emplaceResult.second) {
        CIE_BEGIN_EXCEPTION_TRACING
        this->insert(Vertex(emplaceResult.first.value().source(), {}))
              .mutableEdges().insert(id);
        this->insert(Vertex(emplaceResult.first.value().target(), {}))
              .mutableEdges().insert(id);
        CIE_END_EXCEPTION_TRACING
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
bool Graph<TVD,TED>::eraseVertex(Size id) noexcept
{
    auto it_vertex = _vertices.find(id);
    if (it_vertex != _vertices.end()) {
        Ref<const Vertex> r_vertex = it_vertex.value();

        // Erase associated edges
        for (Size edgeID : r_vertex.edges()) {
            auto it_edge = _edges.find(edgeID);

            // Erase edge from the endpoint vertices' edge lists
            StaticArray<Size,2> endpoints = {it_edge.value().source(), it_edge.value().target()};
            for (Size endpoint : endpoints) {
                if (endpoint != id) {
                    _vertices.find(endpoint).value().mutableEdges().erase(edgeID);
                }
            }

            // Erase edge
            _edges.erase(it_edge);
        }

        // Erase vertex
        _vertices.erase(it_vertex);
        return true;
    } else {
        return false;
    }
}


template <class TVD, class TED>
bool Graph<TVD,TED>::eraseEdge(Size id) noexcept
{
    auto it_edge = _edges.find(id);
    if (it_edge != _edges.end()) {
        Ref<const Edge> r_edge = it_edge.value();

        // Erase edge from the endpoint vertices' edge lists
        _vertices.find(r_edge.source()).value().mutableEdges().erase(id);
        _vertices.find(r_edge.target()).value().mutableEdges().erase(id);

        // Erase edge
        _edges.erase(it_edge);
        return true;
    } else {
        return false;
    }
}


template <class TVD, class TED>
OptionalRef<const typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::findVertex(Size id) const noexcept
{
    return impl::findGraphItem(id, _vertices);
}


template <class TVD, class TED>
OptionalRef<typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::findVertex(Size id) noexcept
{
    return impl::findGraphItem(id, _vertices);
}


template <class TVD, class TED>
OptionalRef<const typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::findEdge(Size id) const noexcept
{
    return impl::findGraphItem(id, _edges);
}


template <class TVD, class TED>
OptionalRef<typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::findEdge(Size id) noexcept
{
    return impl::findGraphItem(id, _edges);
}


} // namespace cie::fem



#endif
