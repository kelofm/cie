#ifndef CIE_FEM_PARTITION_VISITOR_PRE_IMPL_HPP
#define CIE_FEM_PARTITION_VISITOR_PRE_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Sequence.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"


namespace cie::fem::impl {



template <class TAttribute, class ...TAttributes>
void assign(Ref<std::tuple<TAttributes...>> rAttributes, TAttribute&& attribute)
{
    if constexpr (!std::is_same_v<TAttribute,ParentIndex>) {
        std::get<TAttribute>(rAttributes) = std::forward<TAttribute>(attribute);
    }
}


template <std::size_t ContainerIndex, class ...TAttributeContainers>
std::size_t getParentIndex(Ptr<std::size_t> pOutputBegin,
                           std::size_t iAttribute,
                           Ref<const std::tuple<Ptr<const TAttributeContainers>...>> rContainers)
{
    Ref<std::size_t> rParentIndex = pOutputBegin[ContainerIndex];
    rParentIndex = std::get<ContainerIndex>(rContainers)->template at<ParentIndex>(iAttribute);
    return rParentIndex;
}


template <class ...TAttributeContainers, std::size_t ...Is>
void getParentIndices(Ptr<std::size_t> pOutputBegin,
                      std::size_t iLeaf,
                      Ref<const std::tuple<Ptr<const TAttributeContainers>...>> rContainers,
                      ct::IndexSequence<Is...>)
{
    pOutputBegin[ct::PackSize<TAttributeContainers...>] = iLeaf;
    (..., (iLeaf = getParentIndex<Is>(pOutputBegin, iLeaf, rContainers)));
}


template <std::size_t ContainerIndex, class ...TAttributeContainers, class ...TAttributes>
auto getAttributes(Ptr<const std::size_t> pAttributeIndexBegin,
                   Ref<const std::tuple<Ptr<TAttributeContainers>...>> rContainers,
                   Ptr<std::tuple<TAttributes...>>)
{
    return std::tie(
        std::get<ContainerIndex>(rContainers)->template at<TAttributes>(pAttributeIndexBegin[ContainerIndex])...
    );
    //return std::get<ContainerIndex>(rContainers)->get(pAttributeIndexBegin[ContainerIndex]);
}


template <class ...TAttributeContainers, std::size_t ...Is>
auto collectAttributes(Ptr<const std::size_t> pAttributeIndexBegin,
                       Ref<const std::tuple<Ptr<TAttributeContainers>...>> rContainers,
                       ct::IndexSequence<Is...>)
{
    return std::tuple_cat(getAttributes<Is, TAttributeContainers...>(
        pAttributeIndexBegin,
        rContainers,
        Ptr<typename ct::Filter<typename std::tuple_element_t<Is,std::tuple<TAttributeContainers...>>::Values,
                                std::tuple<ParentIndex>>::Type> {}
    )...);
}


struct AttributeAggregate
{
    template <class TAttributeContainer, class ...TAttributeContainers>
    static auto get(Size leafIndex,
                    TAttributeContainer&& rContainer,
                    TAttributeContainers&&... rContainers)
    {
        constexpr auto packSize = ct::PackSize<TAttributeContainers...>;
        std::array<std::size_t,packSize+1> indices;
        const auto immutableContainers = std::make_tuple<Ptr<const std::remove_reference_t<TAttributeContainers>>...>(&rContainers...);
        const auto containers = std::make_tuple(&rContainers...);
        impl::getParentIndices(
            indices.data(),
            leafIndex,
            immutableContainers,
            ct::MakeReverseIndexSequence<packSize>()
        );
        return collectAttributes<
            std::remove_reference_t<TAttributeContainer>,
            std::remove_reference_t<TAttributeContainers>...
        >(indices.data(),
          std::tuple_cat(std::make_tuple(&rContainer), containers),
          ct::MakeIndexSequence<packSize+1>());
    }
}; // struct AttributeAggregate


} // namespace cie::fem::impl


#endif
