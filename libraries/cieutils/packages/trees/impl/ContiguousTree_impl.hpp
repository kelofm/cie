#ifndef CIE_UTILS_CONTIGUOUS_TREE_IMPL_HPP
#define CIE_UTILS_CONTIGUOUS_TREE_IMPL_HPP

// --- Utility Includes ---
#include "packages/trees/inc/ContiguousTree.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <algorithm>
#include <limits>
#include <stack>


namespace cie::utils {


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
ContiguousTree<TNode,TIndex,TTag>::ContiguousTree()
    : _nodes(),
      _availableIndices(),
      _mutex()
{
    _nodes.emplace_back(TNode::makeRoot());
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
void
ContiguousTree<TNode,TIndex,TTag>::erase(TIndex nodeIndex)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(nodeIndex, "Attempt to erase the root node")
    CIE_OUT_OF_RANGE_CHECK(nodeIndex + ContiguousTree::ChildrenPerNode < _nodes.size())
    CIE_DEBUG_CHECK((nodeIndex - 1) % ContiguousTree::ChildrenPerNode == 0, "Attempt to partly erase sibling nodes")

    [[maybe_unused]] const auto lock = _mutex.makeScopedLock();
    std::stack<TIndex,std::vector<TIndex>> indexStack;
    indexStack.push(nodeIndex);

    while (!indexStack.empty()) {
        const auto iNode = indexStack.top();
        indexStack.pop();

        for (unsigned iSibling=0; iSibling<ContiguousTree::ChildrenPerNode; ++iSibling) {
            Ref<TNode> rNode = this->_nodes[iNode + iSibling];
            if (!rNode.isLeaf()) {
                indexStack.push(rNode.getChildBegin());
            }
            rNode.setNull();
        }

        _availableIndices.push(iNode);
    }

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor>
void
ContiguousTree<TNode,TIndex,TTag>::visit(const TFunctor& rVisitFunctor) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(!_nodes.empty())
    using IndexPair = std::pair<TIndex,TIndex>; // {iNode, depth}
    std::stack<IndexPair,DynamicArray<IndexPair>> indexStack;

    if (!_nodes.front().isLeaf()) {
        indexStack.emplace(0, 0);
    }

    while (!indexStack.empty()) {
        const auto indexPair = indexStack.top();
        indexStack.pop();

        for (TIndex iSibling=0; iSibling<TNode::ChildrenPerNode; ++iSibling) {
            const TIndex iNode = indexPair.first + iSibling;
            Ref<const Node> rNode = _nodes[iNode];
            if (!rVisitFunctor(rNode, indexPair.second)) {
                return;
            }
            if (!rNode.isLeaf()) {
                indexStack.emplace(rNode.getChildBegin(), indexPair.second + 1);
            }
        }
    } // while indexStack

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor, concepts::ThreadPool TPool>
void
ContiguousTree<TNode,TIndex,TTag>::visit(const TFunctor&, TPool&) const
{
    CIE_THROW(NotImplementedException, "Multithreaded node visiting is not implemented yet")
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::Iterator<TNode> TIterator>
TIndex ContiguousTree<TNode,TIndex,TTag>::insert(TIterator itNodeBegin)
{
    CIE_BEGIN_EXCEPTION_TRACING

    static_assert(TNode::ChildrenPerNode < std::numeric_limits<TIndex>::max(),
                  "Children cannot be indexed by the specified index type!");

    [[maybe_unused]] const auto lock = _mutex.makeScopedLock();

    TIndex iBegin;

    if (_availableIndices.empty()) {
        // Check for overflows
        const auto size = _nodes.size();
        iBegin = size;
        const auto extendedSize = _nodes.size() + TNode::ChildrenPerNode;
        CIE_CHECK(size < extendedSize && extendedSize <= std::numeric_limits<TIndex>::max(),
                  "Index overflow!")

        _nodes.reserve(extendedSize); // <-- this isn't that relevant here
        for (TIndex i=0; i<ContiguousTree::ChildrenPerNode; ++i, ++itNodeBegin) {
            _nodes.emplace_back(*itNodeBegin);
        }
    } else {
        // Pop an available begin from the set.
        iBegin = _availableIndices.front();
        _availableIndices.pop();

        // Insert nodes
        auto itNode = _nodes.begin() + iBegin;
        const auto itNodeEnd = itNode + ContiguousTree::ChildrenPerNode;
        for (; itNode!=itNodeEnd; ++itNode, ++itNodeBegin) {
            *itNode = *itNodeBegin;
        }
    }

    return iBegin;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
void
ContiguousTree<TNode,TIndex,TTag>::clear() noexcept
{
    auto root = std::move(_nodes[0]);
    root.unsetChildBegin();
    _nodes.clear();
    _nodes.emplace_back(std::move(root));
    _availableIndices = {};
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
void
ContiguousTree<TNode,TIndex,TTag>::reset()
{
    TNode root = std::move(_nodes[0]);
    _nodes = NodeContainer();
    _availableIndices = IndexContainer();
    _nodes.emplace_back(std::move(root));
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
Size
ContiguousTree<TNode,TIndex,TTag>::size() const noexcept
{
    return std::count_if(_nodes.begin(),
                         _nodes.end(),
                         [](const auto& rNode){return !rNode.isNull();});
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
Bool
ContiguousTree<TNode,TIndex,TTag>::empty() const noexcept
{
    return this->size() == 0;
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
Size
ContiguousTree<TNode,TIndex,TTag>::numberOfNodes() const
{
    return _nodes.size();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
Size
ContiguousTree<TNode,TIndex,TTag>::numberOfNulls() const
{
    return _nodes.size() - this->size();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
typename ContiguousTree<TNode,TIndex,TTag>::const_iterator
ContiguousTree<TNode,TIndex,TTag>::begin() const noexcept
{
    return _nodes.begin();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
typename ContiguousTree<TNode,TIndex,TTag>::iterator
ContiguousTree<TNode,TIndex,TTag>::begin() noexcept
{
    return _nodes.begin();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
typename ContiguousTree<TNode,TIndex,TTag>::const_iterator
ContiguousTree<TNode,TIndex,TTag>::end() const noexcept
{
    return _nodes.end();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
typename ContiguousTree<TNode,TIndex,TTag>::iterator
ContiguousTree<TNode,TIndex,TTag>::end() noexcept
{
    return _nodes.end();
}


} // namespace cie::utils


#endif
