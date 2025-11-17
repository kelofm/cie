#ifndef CIE_FEM_PARTITION_VISITOR_HPP
#define CIE_FEM_PARTITION_VISITOR_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/impl/PartitionVisitor_pre_impl.hpp"


namespace cie::fem {


/// @brief Convenience class for fetching attributes from a set of @ref Partition "partitions".
/// @ingroup fem
class PartitionVisitor
{
private:
    template <class ...TPartitions, std::size_t ...Is>
    static auto collectCellAttributes(Size leafIndex,
                                      Ref<const std::tuple<Ptr<TPartitions>...>> rPartitions,
                                      ct::IndexSequence<Is...>)
    {
        return impl::AttributeAggregate::get(
            leafIndex,
            std::get<Is>(rPartitions)->_cellAttributes ...
        );
    }

public:
    /// @brief  Fetch all attributes of an entry defined by its index in the leaf partition.
    /// @tparam TPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @param leafIndex The entry's index in the leaf partition.
    /// @param ...rPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @return A tuple with all attributes, ordered from from root to leaf (top to bottom).
    template <class ...TPartitions>
    requires (concepts::Partition<TPartitions> && ...)
    static auto collectCellAttributes(Size leafIndex, Ref<const TPartitions>... rPartitions)
    {
        return impl::AttributeAggregate::get(
            leafIndex,
            rPartitions._cellAttributes ...
        );
    }


    /// @brief  Fetch all attributes of an entry defined by its index in the leaf partition.
    /// @tparam TPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @param leafIndex The entry's index in the leaf partition.
    /// @param ...pPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @return A tuple with all attributes, ordered from from root to leaf (top to bottom).
    template <class ...TPartitions>
    requires (concepts::Partition<TPartitions> && ...)
    static auto collectCellAttributes(Size leafIndex, const TPartitions* ... pPartitions)
    {
        return PartitionVisitor::collectCellAttributes<TPartitions...>(leafIndex, *pPartitions ...);
    }


    /// @brief  Fetch all attributes of an entry defined by its index in the leaf partition.
    /// @tparam TPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @param leafIndex The entry's index in the leaf partition.
    /// @param ...rPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @return A tuple with all attributes, ordered from from root to leaf (top to bottom).
    template <class ...TPartitions>
    requires (concepts::Partition<TPartitions> && ...)
    static auto collectCellAttributes(Size leafIndex, Ref<const std::tuple<Ptr<const TPartitions>...>> rPartitions)
    {
        return PartitionVisitor::collectCellAttributes(leafIndex, rPartitions, ct::MakeIndexSequence<ct::PackSize<TPartitions...>>());
    }

    /// @brief  Fetch all attributes of an entry defined by its index in the leaf partition.
    /// @tparam ...TPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @param leafIndex The entry's index in the leaf partition.
    /// @param ...rPartitions Set of partitions to fetch attributes from. First is root, last is leaf,
    ///         and every partition between them must be ordered accordingly.
    /// @return A tuple with all attributes, ordered from from root to leaf (top to bottom).
    template <class ...TPartitions>
    requires (concepts::Partition<TPartitions> && ...)
    static auto collectCellAttributes(Size leafIndex, Ref<const std::tuple<Ptr<TPartitions>...>> rPartitions)
    {
        return PartitionVisitor::collectCellAttributes(
            leafIndex,
            rPartitions,
            ct::MakeIndexSequence<ct::PackSize<TPartitions...>>());
    }
}; // class PartitionVisitor


} // namespace cie::fem


#endif
