#ifndef CIE_FEM_ASSEMBLER_HPP
#define CIE_FEM_ASSEMBLER_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"

// --- STL Includes ---
#include <packages/macros/inc/checks.hpp>
#include <unordered_map>
#include <vector>
#include <optional>
#include <ranges> // transform_view


namespace cie::fem {


class Assembler
{
private:
    using DefaultMesh = Graph<void,void>;

    using DoFMap = std::unordered_map<
        DefaultMesh::VertexID,                  // <== ID of the cell
        std::vector<std::optional<std::size_t>> // <== DoF indices of its basis functions
    >;

    static auto makeIndexView(Ref<const DoFMap::mapped_type> rIndices)
    {
        return std::ranges::transform_view(
            rIndices,
            [](DoFMap::mapped_type::value_type maybeIndex){return maybeIndex.value();}
        );
    }

public:
    using DoFPairVector = std::vector<std::pair<std::size_t,std::size_t>>;

    using DoFPairIterator = std::back_insert_iterator<DoFPairVector>;

    Assembler() noexcept;

    explicit Assembler(std::size_t dofBegin) noexcept;

    template <class TVertexData,
              class TEdgeData,
              concepts::FunctionWithSignature<std::size_t,Ref<const typename Graph<TVertexData,TEdgeData>::Vertex>> TDoFCounter,
              concepts::FunctionWithSignature<void,Ref<const typename Graph<TVertexData,TEdgeData>::Edge>,DoFPairIterator> TDoFPairFunctor>
    void addGraph(Ref<const Graph<TVertexData,TEdgeData>> rGraph,
                  TDoFCounter&& rDoFCounter,
                  TDoFPairFunctor&& rDoFMatcher);

    std::size_t dofCount() const noexcept;

    auto keys() const
    {
        return std::ranges::transform_view(
            _dofMap,
            [](const auto& rPair) {return rPair.first;}
        );
    }

    auto items() const
    {
        return std::ranges::transform_view(
            _dofMap,
            [](const auto& rPair) {
                return std::make_pair(rPair.first,
                                      Assembler::makeIndexView(rPair.second));
            }
        );
    }

    auto operator[](DefaultMesh::VertexID vertexID) const
    {
        const auto it = _dofMap.find(vertexID);
        CIE_OUT_OF_RANGE_CHECK(it != _dofMap.end())
        return this->makeIndexView(it->second);
    }

private:
    std::size_t _dofCounter;

    DoFMap _dofMap;
}; // class Assembler


} // namespace cie::fem

#include "packages/graph/impl/Assembler_impl.hpp"

#endif
