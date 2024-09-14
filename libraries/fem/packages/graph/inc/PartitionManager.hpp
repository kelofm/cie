#ifndef CIE_FEM_PARTITION_MANAGER_HPP
#define CIE_FEM_PARTITION_MANAGER_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/impl/PartitionManager_pre_impl.hpp"


namespace cie::fem {


class PartitionManager
{
public:
    template <PartitionBase::Item TItem, class ...TPartitions>
    static auto collectAttributes(Size leafIndex, Ref<const TPartitions>... rPartitions)
    {
        if constexpr (TItem == PartitionBase::Item::Vertex) {
            return impl::AttributeAggregate<typename TPartitions::VertexAttributes...>::get(
                rPartitions._vertexAttributes ...,
                leafIndex
            );
        } else if constexpr (TItem == PartitionBase::Item::Polytope) {
            return impl::AttributeAggregate<typename TPartitions::PolytopeAttributes...>::get(
                rPartitions._polytopeAttributes ...,
                leafIndex
            );
        } else {
            static_assert(std::is_same_v<std::tuple<TPartitions...>,std::tuple<>>);
        }
    }
}; // class PartitionManager


} // namespace cie::fem


#endif
