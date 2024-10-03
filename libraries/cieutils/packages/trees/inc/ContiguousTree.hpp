#ifndef CIE_UTILS_CONTIGUOUS_TREE_HPP
#define CIE_UTILS_CONTIGUOUS_TREE_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
#include "packages/types/inc/tags.hpp"
#include "packages/concurrency/inc/Mutex.hpp"

// --- STL Includes ---
#include <iterator>
#include <queue>


namespace cie::concepts {


/// @addtogroup cieutils
template <class T>
concept ContiguousTreeNode
= requires (T instance, const T constInstance)
{
    typename T::Index; // index type
    T::ChildrenPerNode; // unsigned integer value

    {T::makeRoot()} -> std::same_as<T>;
    {constInstance.isLeaf()} -> std::same_as<Bool>;
    {constInstance.getChildBegin()} -> std::same_as<typename T::Index>;
    {constInstance.getParentIndex()} -> std::same_as<typename T::Index>;
    {instance.unsetChildBegin()};
    {instance.setNull()};
    {constInstance.isNull()} -> std::same_as<Bool>;
}; // concept ContiguousTreeNode


} // namespace cie::concepts


namespace cie::utils {


/// @addtogroup cieutils
template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TParallelTag = tags::Serial>
class ContiguousTree
{
private:
    using NodeContainer = DynamicArray<TNode>;

    using IndexContainer = std::queue<TIndex>;

    using Mutex = mp::Mutex<TParallelTag>;

    static inline constexpr const Size ChildrenPerNode = TNode::ChildrenPerNode;

public:
    CIE_DEFINE_CLASS_POINTERS(ContiguousTree)

    using Node = TNode;

    using Index = TIndex;

    using iterator = typename NodeContainer::iterator;

    using const_iterator = typename NodeContainer::const_iterator;

    using value_type = TNode;

    using size_type = typename NodeContainer::size_type;

public:
    ContiguousTree();

    CIE_DEFINE_CLASS_DEFAULT_MOVES(ContiguousTree)

    CIE_DEFINE_CLASS_DEFAULT_COPIES(ContiguousTree)

    /// @brief Erase @ref TNode::ChildrenPerNode nodes at once (thread safe if tags::SMP) beginning at @a nodeIndex.
    void erase(TIndex nodeIndex);

    /// @brief Invoke a functor on every node in the tree.
    /// @details The callable functor should take a reference to a node and
    ///          an index representing its depth. If the depth is not required,
    ///          use @ref begin and @ref end to loop through the nodes in
    ///          unspecified order instead.
    /// @param rVisitFunctor: function to invoke on each node. It must take a node
    ///                        and level as argument, and produce a bool that
    ///                        indicates whether the visit loop should be continued.
    /// @details Single thread version.
    template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor>
    void visit(const TFunctor& rVisitFunctor) const;

    /// @brief Invoke a functor on every node in the tree.
    /// @details The callable functor should take a reference to a node and
    ///          an index representing its depth. If the depth is not required,
    ///          use @ref begin and @ref end to loop through the nodes in
    ///          unspecified order instead.
    /// @param rVisitFunctor: function to invoke on each node. It must take a node
    ///                        and level as argument, and produce a bool that
    ///                        indicates whether the visit loop should be continued.
    /// @param rPool: Thread pool to use for multithreading.
    /// @details Multithread version.
    /// @todo Implement multithreaded version.
    template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor, concepts::ThreadPool TPool>
    void visit(const TFunctor& rVisitFunctor, TPool& rPool) const;

    /// @brief Insert @ref TNode::ChildrenPerNode nodes at once (thread safe if tags::SMP).
    /// @return Index of the first inserted node in the internal node list.
    template <concepts::Iterator<TNode> TIterator>
    TIndex insert(TIterator itNodeBegin);

    /// @brief Reset the tree (all iterators and indices are invalidated).
    void clear() noexcept;

    /// @brief Free allocated resources but keep the root node.
    void reset();

    Size size() const noexcept;

    Bool empty() const noexcept;

    const_iterator begin() const noexcept;

    iterator begin() noexcept;

    const_iterator end() const noexcept;

    iterator end() noexcept;

private:
    /// @brief Get the total number of nodes, including nulls.
    Size numberOfNodes() const;

    Size numberOfNulls() const;

protected:
    NodeContainer _nodes;

private:
    /// @brief Indices pointing to the beginning of @ref TNode::ChildrenPerNode long sets of erased nodes.
    IndexContainer _availableIndices;

    mutable Mutex _mutex;
}; // class ContiguousTree


} // cie::utils

#include "packages/trees/impl/ContiguousTree_impl.hpp"

#endif
