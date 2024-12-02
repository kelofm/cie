#ifndef CIE_FEM_PARTITION_IMPL_HPP
#define CIE_FEM_PARTITION_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"


namespace cie::fem {


template <concepts::AttributeContainer TCAC,
          class TPA>
Partition<TCAC,TPA>::Partition(RightRef<std::string> rName,
                               PartitionID id)
    : Partition(std::move(rName),
                id,
                {})
{
}


template <concepts::AttributeContainer TCAC,
          class TPA>
Partition<TCAC,TPA>::Partition(RightRef<std::string> rName,
                               PartitionID id,
                               RightRef<CellAttributeContainer> rCellAttributes) noexcept
    requires std::is_same_v<Attributes,std::monostate>
    : PartitionBase(std::move(rName), id),
      _cellAttributes(std::move(rCellAttributes)),
      _attributes()
{
}


template <concepts::AttributeContainer TCAC,
          class TPA>
Partition<TCAC,TPA>::Partition(RightRef<std::string> rName,
                               PartitionID id,
                               RightRef<CellAttributeContainer> rCellAttributes,
                               RightRef<Attributes> rPartitionAttributes) noexcept
    requires (!std::is_same_v<Attributes,std::monostate>)
    : PartitionBase(std::move(rName), id),
      _cellAttributes(std::move(rCellAttributes)),
      _attributes(std::move(rPartitionAttributes))
{
}


template <concepts::Partition TPartition, class ...TChildren>
template <unsigned ...Is>
Ref<const typename PartitionTree<TPartition,TChildren...>::template Child<Is...>>
PartitionTree<TPartition,TChildren...>::child() const
{
    if constexpr (sizeof ...(Is) != 0) {
        return *ChildSelector<Is...>::get(this);
    } else {
        return *this;
    }
}


template <concepts::Partition TPartition, class ...TChildren>
template <unsigned ...Is>
Ref<typename PartitionTree<TPartition,TChildren...>::template Child<Is...>>
PartitionTree<TPartition,TChildren...>::child()
{
    if constexpr (sizeof ...(Is) != 0) {
        return *ChildSelector<Is...>::get(this);
    } else {
        return *this;
    }
}


template <concepts::Partition TPartition, class ...TChildren>
template <unsigned ...Is>
requires (sizeof...(Is) != 0)
Ref<typename PartitionTree<TPartition,TChildren...>::template Child<Is...>>
PartitionTree<TPartition,TChildren...>::setChild(RightRef<std::unique_ptr<Child<Is...>>> rpChild)
{
    return *(ChildSelector<Is...>::get(this) = std::move(rpChild));
}


template <unsigned ...NestedIndices, concepts::Partition TRoot>
PartitionHandle<TRoot,NestedIndices...> makePartitionHandle(Ref<TRoot> rRoot) noexcept
{
    return PartitionHandle<TRoot,NestedIndices...>(rRoot);
}


template <unsigned ...NestedIndices, concepts::Partition TRoot>
PartitionHandle<TRoot,NestedIndices...> makePartitionHandle(Ref<TRoot> rRoot,
                                                            std::integer_sequence<unsigned,NestedIndices...>) noexcept
{
    return PartitionHandle<TRoot,NestedIndices...>(rRoot);
}


template <concepts::Partition TRoot, unsigned ...Is>
PartitionHandle<TRoot,Is...>::PartitionHandle(Ref<TRoot> rRoot)
    : _pRoot(&rRoot),
      _pPartition(&rRoot.template child<Is...>())
{
}


template <concepts::Partition TRoot, unsigned ...NestedIndices>
Ref<const typename PartitionHandle<TRoot,NestedIndices...>::element_type>
PartitionHandle<TRoot,NestedIndices...>::operator*() const noexcept
{
    return *_pPartition;
}


template <concepts::Partition TRoot, unsigned ...NestedIndices>
Ref<typename PartitionHandle<TRoot,NestedIndices...>::element_type>
PartitionHandle<TRoot,NestedIndices...>::operator*() noexcept
{
    return *_pPartition;
}


template <concepts::Partition TRoot, unsigned ...NestedIndices>
Ref<const TRoot>
PartitionHandle<TRoot,NestedIndices...>::root() const noexcept
{
    return *_pPartition;
}


template <concepts::Partition TRoot, unsigned ...NestedIndices>
Ref<TRoot>
PartitionHandle<TRoot,NestedIndices...>::root() noexcept
{
    return *_pPartition;
}


} // namespace cie::fem


#endif
