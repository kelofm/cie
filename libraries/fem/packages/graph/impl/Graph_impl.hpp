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


template <class TID, class TContainer>
OptionalRef<typename CopyConstQualifier<
    TContainer,
    typename std::remove_reference_t<TContainer>::mapped_type
>::Type> findGraphItem(TID id, TContainer&& r_container) noexcept
{
    const auto it = r_container.find(id);
    if (it == r_container.end()) {
        return {};
    } else {
        return it.value();
    }
}


template <class TSecond, class TFirst>
std::conditional_t<
    std::is_same_v<std::remove_const_t<TSecond>,void>,
    std::tuple<TFirst>,
    std::tuple<TFirst,TSecond>
> makePartiallyInitializedTuple(TFirst&& r_first)
{
    if constexpr (std::is_same_v<std::remove_const_t<TSecond>,void>) {
        return std::tuple<TFirst>(std::forward<TFirst>(r_first));
    } else {
        return std::tuple<TFirst,TSecond>(std::forward<TFirst>(r_first), TSecond());
    }
}


} // namespace impl


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(VertexID id,
                               RightRef<tsl::robin_set<EdgeID>> rEdges) noexcept
    : ItemBase<VertexID>(id),
      _data(impl::makePartiallyInitializedTuple<TVD>(std::move(rEdges)))
{}


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(VertexID id,
                               RightRef<tsl::robin_set<EdgeID>> rEdges,
                               std::conditional_t<
                                    std::is_same_v<std::remove_const_t<TVD>,void>,
                                    char, // dummy type, cannot be void
                                    typename VoidSafe<TVD>::RightRef
                               > rData) noexcept
requires (!std::is_same_v<std::remove_const_t<TVD>,void>)
    : ItemBase<VertexID>(id),
      _data(std::move(rEdges), std::move(rData))
{}


template <class TVD, class TED>
Graph<TVD,TED>::Vertex::Vertex(VertexID id,
                               RightRef<tsl::robin_set<EdgeID>> rEdges,
                               std::conditional_t<
                                    std::is_same_v<std::remove_const_t<TVD>,void>,
                                    int, // dummy type, cannot be void
                                    typename VoidSafe<const TVD>::Ref
                               > rData)
requires (!std::is_same_v<std::remove_const_t<TVD>,void>)
    : ItemBase<VertexID>(id),
      _data(std::move(rEdges), rData)
{}


template <class TVD, class TED>
Ref<const tsl::robin_set<typename Graph<TVD,TED>::EdgeID>>
Graph<TVD,TED>::Vertex::edges() const noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
typename VoidSafe<const TVD>::Ref Graph<TVD,TED>::Vertex::data() const noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TVD>,void>) {
        return std::get<1>(_data);
    }
}


template <class TVD, class TED>
typename VoidSafe<TVD>::Ref Graph<TVD,TED>::Vertex::data() noexcept
{
    if constexpr (!std::is_same_v<std::remove_const_t<TVD>,void>) {
        return std::get<1>(_data);
    }
}


template <class TVD, class TED>
Ref<tsl::robin_set<typename Graph<TVD,TED>::EdgeID>> Graph<TVD,TED>::Vertex::mutableEdges() noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(EdgeID id,
                           std::pair<VertexID,VertexID> vertices) noexcept
    : ItemBase<EdgeID>(id),
      _data(impl::makePartiallyInitializedTuple<TED>(vertices))
{
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(EdgeID id,
                           std::pair<VertexID,VertexID> vertices,
                           std::conditional_t<
                                std::is_same_v<std::remove_const_t<TED>,void>,
                                char, // dummy type, cannot be void
                                typename VoidSafe<TED>::RightRef
                           > rData) noexcept
requires (!std::is_same_v<std::remove_const_t<TED>,void>)
    : ItemBase<EdgeID>(id),
      _data(vertices, std::move(rData))
{
}


template <class TVD, class TED>
Graph<TVD,TED>::Edge::Edge(EdgeID id,
                           std::pair<VertexID,VertexID> vertices,
                           std::conditional_t<
                                std::is_same_v<std::remove_const_t<TED>,void>,
                                int, // dummy type, cannot be void
                                typename VoidSafe<const TED>::Ref
                           > rData)
requires (!std::is_same_v<std::remove_const_t<TED>,void>)
    : ItemBase<EdgeID>(id),
      _data(vertices, rData)
{
}


template <class TVD, class TED>
std::pair<typename Graph<TVD,TED>::VertexID,typename Graph<TVD,TED>::VertexID>
Graph<TVD,TED>::Edge::vertices() const noexcept
{
    return std::get<0>(_data);
}


template <class TVD, class TED>
typename Graph<TVD,TED>::VertexID
Graph<TVD,TED>::Edge::source() const noexcept
{
    return std::get<0>(_data).first;
}


template <class TVD, class TED>
typename Graph<TVD,TED>::VertexID
Graph<TVD,TED>::Edge::target() const noexcept
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
Ref<typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::insert(RightRef<Vertex> rVertex,
                                                            bool overwrite)
{
    CIE_CHECK(rVertex.edges().empty(), "Attempt to insert vertex " << rVertex.id() << " that already has edges")

    const auto id = rVertex.id();
    std::pair<typename decltype(_vertices)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no vertex exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseVertex(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _vertices.emplace(id, std::move(rVertex));
    CIE_END_EXCEPTION_TRACING

    // Check whether every edge connected to
    // the newly inserted vertex exists.
    if (emplaceResult.second) {
        for (EdgeID edgeID : emplaceResult.first.value().edges()) {
            CIE_CHECK(
                this->findEdge(edgeID).has_value(),
                "Vertex " << id << " is connected to edge " << edgeID << " that is not part of the graph\n"
            );
        }
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
Ref<typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::insert(Ref<const Vertex> rVertex,
                                                            bool overwrite)
{
    CIE_CHECK(rVertex.edges().empty(), "Attempt to insert vertex " << rVertex.id() << " that already has edges")

    const auto id = rVertex.id();
    std::pair<typename decltype(_vertices)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no vertex exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseVertex(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _vertices.emplace(id, rVertex);
    CIE_END_EXCEPTION_TRACING

    // Check whether every edge connected to
    // the newly inserted vertex exists.
    if (emplaceResult.second) {
        for (EdgeID edgeID : emplaceResult.first.value().edges()) {
            CIE_CHECK(
                this->findEdge(edgeID).has_value(),
                "Vertex " << id << " is connected to edge " << edgeID << " that is not part of the graph\n"
            );
        }
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
Ref<typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::insert(RightRef<Edge> rEdge,
                                                          bool overwrite)
{
    const auto id = rEdge.id();
    std::pair<typename decltype(_edges)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no edge exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseEdge(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _edges.emplace(id, std::move(rEdge));
    CIE_END_EXCEPTION_TRACING

    // Ensure that both endpoints of the edge exist in the graph,
    // and both end points reference the edge.
    if (emplaceResult.second) {
        CIE_BEGIN_EXCEPTION_TRACING
        this->insert(Vertex(emplaceResult.first.value().source(), {}))
              .mutableEdges()
              .insert(id);
        this->insert(Vertex(emplaceResult.first.value().target(), {}))
              .mutableEdges()
              .insert(id);
        CIE_END_EXCEPTION_TRACING
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
Ref<typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::insert(Ref<const Edge> rEdge,
                                                          bool overwrite)
{
    const auto id = rEdge.id();
    std::pair<typename decltype(_edges)::iterator,bool> emplaceResult {{}, false};

    // Make sure that no edge exists in the graph with the given ID
    // if overwriting was requested.
    if (overwrite) {
        this->eraseEdge(id);
    }

    CIE_BEGIN_EXCEPTION_TRACING
    emplaceResult = _edges.emplace(id, rEdge);
    CIE_END_EXCEPTION_TRACING

    // Ensure that both endpoints of the edge exist in the graph,
    // and both end points reference the edge.
    if (emplaceResult.second) {
        CIE_BEGIN_EXCEPTION_TRACING
        this->insert(Vertex(emplaceResult.first.value().source(), {}))
              .mutableEdges()
              .insert(id);
        this->insert(Vertex(emplaceResult.first.value().target(), {}))
              .mutableEdges()
              .insert(id);
        CIE_END_EXCEPTION_TRACING
    }

    return emplaceResult.first.value();
}


template <class TVD, class TED>
bool Graph<TVD,TED>::eraseVertex(VertexID id) noexcept
{
    auto it_vertex = _vertices.find(id);
    if (it_vertex != _vertices.end()) {
        Ref<const Vertex> rVertex = it_vertex.value();

        // Erase associated edges
        for (EdgeID edgeID : rVertex.edges()) {
            auto it_edge = _edges.find(edgeID);

            // Erase edge from the endpoint vertices' edge lists
            StaticArray<VertexID,2> endpoints = {it_edge.value().source(), it_edge.value().target()};
            for (VertexID endpoint : endpoints) {
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
bool Graph<TVD,TED>::eraseEdge(EdgeID id) noexcept
{
    auto it_edge = _edges.find(id);
    if (it_edge != _edges.end()) {
        Ref<const Edge> rEdge = it_edge.value();

        // Erase edge from the endpoint vertices' edge lists
        _vertices.find(rEdge.source()).value().mutableEdges().erase(id);
        _vertices.find(rEdge.target()).value().mutableEdges().erase(id);

        // Erase edge
        _edges.erase(it_edge);
        return true;
    } else {
        return false;
    }
}


template <class TVD, class TED>
OptionalRef<const typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::findVertex(VertexID id) const noexcept
{
    return impl::findGraphItem(id, _vertices);
}


template <class TVD, class TED>
OptionalRef<typename Graph<TVD,TED>::Vertex> Graph<TVD,TED>::findVertex(VertexID id) noexcept
{
    return impl::findGraphItem(id, _vertices);
}


template <class TVD, class TED>
OptionalRef<const typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::findEdge(EdgeID id) const noexcept
{
    return impl::findGraphItem(id, _edges);
}


template <class TVD, class TED>
OptionalRef<typename Graph<TVD,TED>::Edge> Graph<TVD,TED>::findEdge(EdgeID id) noexcept
{
    return impl::findGraphItem(id, _edges);
}


} // namespace cie::fem



#endif
