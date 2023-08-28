#ifndef CIE_FEM_PARTITION_MANAGER_PRE_IMPL_HPP
#define CIE_FEM_PARTITION_MANAGER_PRE_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Sequence.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"


namespace cie::fem::impl {



template <class TAttribute, class ...TAttributes>
void assign(Ref<std::tuple<TAttributes...>> r_attributes, TAttribute&& attribute)
{
    if constexpr (!std::is_same_v<TAttribute,ParentIndex>) {
        std::get<TAttribute>(r_attributes) = std::forward<TAttribute>(attribute);
    }
}


template <std::size_t ContainerIndex, class ...TAttributeContainers>
std::size_t getParentIndex(Ptr<std::size_t> p_outputBegin,
                           std::size_t i_attribute,
                           Ref<const std::tuple<Ptr<const TAttributeContainers>...>> r_containers)
{
    Ref<std::size_t> r_parentIndex = *(p_outputBegin + ContainerIndex);
    r_parentIndex = std::get<ContainerIndex>(r_containers)->template at<ParentIndex>(i_attribute);
    return r_parentIndex;
}


template <class ...TAttributeContainers, std::size_t ...Is>
void getParentIndices(Ptr<std::size_t> p_outputBegin,
                      std::size_t i_leaf,
                      Ref<const std::tuple<Ptr<const TAttributeContainers>...>> r_containers,
                      ct::IndexSequence<Is...>)
{
    *(p_outputBegin + ct::PackSize<TAttributeContainers...>) = i_leaf;
    (..., (i_leaf = getParentIndex<Is>(p_outputBegin, i_leaf, r_containers)));
}


template <std::size_t ContainerIndex, class ...TAttributeContainers, class ...TAttributes>
auto getAttributes(Ptr<const std::size_t> p_attributeIndexBegin,
                   Ref<const std::tuple<Ptr<const TAttributeContainers>...>> r_containers,
                   Ptr<std::tuple<TAttributes...>>)
{
    return std::make_tuple(
        std::get<ContainerIndex>(r_containers)->template at<TAttributes>(p_attributeIndexBegin[ContainerIndex])...
    );
    //return std::get<ContainerIndex>(r_containers)->get(p_attributeIndexBegin[ContainerIndex]);
}


template <class ...TAttributeContainers, std::size_t ...Is>
auto collectAttributes(Ptr<const std::size_t> p_attributeIndexBegin,
                       Ref<const std::tuple<Ptr<const TAttributeContainers>...>> r_containers,
                       ct::IndexSequence<Is...>)
{
    return std::tuple_cat(getAttributes<Is>(
        p_attributeIndexBegin,
        r_containers,
        Ptr<typename ct::Filter<typename std::tuple_element_t<Is,std::tuple<TAttributeContainers...>>::Values,
                                std::tuple<ParentIndex>>::Type> {}
    )...);
}


template <class TAttributeContainer, class ...TAttributeContainers>
struct AttributeAggregate
{
    static auto get(const TAttributeContainer& r_container,
                    const TAttributeContainers&... r_containers,
                    Size leafIndex)
    {
        constexpr auto packSize = ct::PackSize<TAttributeContainers...>;
        std::array<std::size_t,packSize+1> indices;
        const auto containers = std::make_tuple(&r_containers...);
        impl::getParentIndices(
            indices.data(),
            leafIndex,
            containers,
            ct::MakeReverseIndexSequence<packSize>()
        );
        return collectAttributes(indices.data(),
                                 std::tuple_cat(std::make_tuple(&r_container), containers),
                                 ct::MakeIndexSequence<packSize+1>());
    }
}; // struct AttributeAggregate



} // namespace cie::fem::impl


#endif
