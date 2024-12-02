#ifndef CIE_FEM_PARTITION_HPP
#define CIE_FEM_PARTITION_HPP

// --- FEM Includes ---
#include "packages/utilities/inc/AttributeContainer.hpp"

// --- Utility Includes ---
#include "packages/types/inc/NamedObject.hpp"
#include "packages/types/inc/IDObject.hpp"
#include "packages/stl_extension/inc/StrongTypeDef.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Select.hpp"

// --- STL Includes ---
#include <concepts>
#include <memory> // std::shared_ptr, std::weak_ptr
#include <variant> // std::monostate
#include <tuple> // std::tuple


namespace cie::fem {
class PartitionBase;
} // namespace cie::fem


namespace cie::concepts {


/// @ingroup fem
template <class T>
concept Partition
= std::derived_from<T,fem::PartitionBase> && requires () {
    typename T::Attributes;
    typename T::CellAttributeContainer;
};


/// @ingroup fem
template <class T>
concept MaybePartition
= Partition<T> || std::is_same_v<T,std::monostate>;


} // namespace cie::concepts


namespace cie::fem {


/// @ingroup fem
CIE_STRONG_TYPEDEF(unsigned, PartitionID);


/// @ingroup fem
CIE_STRONG_TYPEDEF(Size, ParentIndex);


/// @ingroup fem
class PartitionBase : public utils::NamedObject,
                      public utils::IDObject<PartitionID>
{
public:
    using ID = PartitionID;

public:
    PartitionBase() noexcept;

    PartitionBase(const PartitionBase&) = delete;

    PartitionBase(PartitionBase&&) noexcept = default;

    PartitionBase(RightRef<std::string> rName, PartitionID id);

    PartitionBase& operator=(const PartitionBase&) = delete;

    PartitionBase& operator=(PartitionBase&&) noexcept = default;

    virtual ~PartitionBase() = default;
}; // class PartitionBase



/// @ingroup fem
template <concepts::AttributeContainer TCellAttributeContainer,
          class TPartitionAttributes = std::monostate>
class Partition : public PartitionBase
{
public:
    using Attributes = TPartitionAttributes;

    using CellAttributeContainer = TCellAttributeContainer;

public:
    Partition() noexcept = default;

    Partition(Partition&& rRhs) noexcept = default;

    Partition(const Partition& rRhs) = delete;

    Partition(RightRef<std::string> rName,
              PartitionID id);

    Partition(RightRef<std::string> rName,
              PartitionID id,
              RightRef<CellAttributeContainer> rCellAttributes) noexcept
        requires std::is_same_v<Attributes,std::monostate>;

    Partition(RightRef<std::string> rName,
              PartitionID id,
              RightRef<CellAttributeContainer> rCellAttributes,
              RightRef<Attributes> rPartitionAttributes) noexcept
        requires (!std::is_same_v<Attributes,std::monostate>);

    Partition& operator=(Partition&& rRhs) noexcept = default;

    Partition& operator=(const Partition& rRhs) = delete;

    Ref<const Attributes> attributes() const noexcept
    {return _attributes;}

    Ref<Attributes> attributes() noexcept
    {return _attributes;}

private:
    friend class PartitionVisitor;

    TCellAttributeContainer _cellAttributes;

    Attributes _attributes;
}; // class Partition


/// @ingroup fem
template <concepts::Partition TPartition, class ...TChildren>
class PartitionTree : public TPartition
{
private:
    template <concepts::Partition TP, class ...TC>
    friend class PartitionTree;

    template <unsigned ...Indices>
    struct ChildSelector {using Type = void;};

    template <unsigned I> struct ChildSelector<I>
    {
        using Type = typename ct::Select<TChildren...>:: template At<I>;
        static Ref<const std::unique_ptr<Type>> get(Ptr<const PartitionTree> pPartition) {
            return std::get<I>(pPartition->_children);
        }
        static Ref<std::unique_ptr<Type>> get(Ptr<PartitionTree> pPartition) {
            return std::get<I>(pPartition->_children);
        }
    };

    template <unsigned I, unsigned ...Is> struct ChildSelector<I,Is...>
    {
        using Type = typename ct::Select<TChildren...>::template At<I>::template ChildSelector<Is...>::Type;
        static Ref<const std::unique_ptr<Type>> get(Ptr<const PartitionTree> pPartition) {
            using SubSelector = typename std::tuple_element_t<I,std::tuple<TChildren...>>::template ChildSelector<Is...>;
            return SubSelector::get(&pPartition->child<I>());
        }
        static Ref<std::unique_ptr<Type>> get(Ptr<PartitionTree> pPartition) {
            using SubSelector = typename std::tuple_element_t<I,std::tuple<TChildren...>>::template ChildSelector<Is...>;
            return SubSelector::get(&pPartition->child<I>());
        }
    };

public:
    using Base = TPartition;

    constexpr static unsigned ChildCount = sizeof...(TChildren);

    using Children = std::tuple<TChildren...>;

    using ChildContainer = std::tuple<std::unique_ptr<TChildren>...>;

    template <unsigned ...Is>
    using Child = std::conditional_t<
        sizeof ...(Is) != 0,
        typename ChildSelector<Is...>::Type,
        TPartition
    >;

public:
    PartitionTree() noexcept = default;

    PartitionTree(PartitionTree&& rRhs) noexcept = default;

    PartitionTree(const PartitionTree& rRhs) = delete;

    using TPartition::TPartition;

    PartitionTree& operator=(PartitionTree&& rRhs) noexcept = default;

    PartitionTree& operator=(const PartitionTree& rRhs) = delete;

    template <unsigned ...Is>
    Ref<const Child<Is...>> child() const;

    template <unsigned ...Is>
    Ref<Child<Is...>> child();

    template <unsigned ...Is>
    requires (sizeof...(Is) != 0)
    Ref<Child<Is...>> setChild(RightRef<std::unique_ptr<Child<Is...>>> rpChild);

private:
    friend class PartitionVisitor;

    ChildContainer _children;
}; // class PartitionTree


template <concepts::Partition TRoot, unsigned ...NestedIndices>
class PartitionHandle;


/// @ingroup fem
template <unsigned ...NestedIndices, concepts::Partition TRoot>
PartitionHandle<TRoot,NestedIndices...> makePartitionHandle(Ref<TRoot> rRoot) noexcept;


/// @ingroup fem
template <unsigned ...NestedIndices, concepts::Partition TRoot>
PartitionHandle<TRoot,NestedIndices...> makePartitionHandle(Ref<TRoot> rRoot,
                                                            std::integer_sequence<unsigned,NestedIndices...>) noexcept;


/// @ingroup fem
template <concepts::Partition TRoot, unsigned ...NestedIndices>
class PartitionHandle
{
public:
    using pointer = PartitionHandle;

    using element_type = typename TRoot::template Child<NestedIndices...>;

    using difference_type = void;

public:
    PartitionHandle() noexcept = default;

    explicit PartitionHandle(Ref<TRoot> rRoot);

    Ref<const element_type> operator*() const noexcept;

    Ref<element_type> operator*() noexcept;

    Ref<const TRoot> root() const noexcept;

    Ref<TRoot> root() noexcept;

    auto flatten() const
    requires (sizeof...(NestedIndices) != 0)
    {
        using IntegralSelector = ct::SelectIntegral<unsigned,NestedIndices...>;
        return std::tuple_cat(
            std::make_tuple(_pRoot),
            makePartitionHandle(_pRoot->template child<IntegralSelector::First>(),
                                typename IntegralSelector::Bottom()).flatten()
        );
    }

    auto flatten()
    requires (sizeof...(NestedIndices) != 0)
    {
        using IntegralSelector = ct::SelectIntegral<unsigned,NestedIndices...>;
        return std::tuple_cat(
            std::make_tuple(_pRoot),
            makePartitionHandle(_pRoot->template child<IntegralSelector::First>(),
                                typename IntegralSelector::Bottom()).flatten()
        );
    }

    auto flatten() const
    requires (sizeof...(NestedIndices) == 0)
    {
        return std::make_tuple(_pRoot);
    }

    auto flatten()
    requires (sizeof...(NestedIndices) == 0)
    {
        return std::make_tuple(_pRoot);
    }

private:
    friend class PartitionVisitor;

    Ptr<TRoot> _pRoot;

    Ptr<element_type> _pPartition;
}; // class PartitionHandle


} // namespace cie::fem

#include "packages/graph/impl/Partition_impl.hpp"

#endif
