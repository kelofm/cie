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
inline void
ContiguousTree<TNode,TIndex,TTag>::erase(TIndex nodeIndex)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(nodeIndex, "Attempt to erase the root node")
    CIE_OUT_OF_RANGE_CHECK(nodeIndex + ContiguousTree::ChildrenPerNode < _nodes.size())
    CIE_DEBUG_CHECK((nodeIndex - 1) % ContiguousTree::ChildrenPerNode == 0, "Attempt to partly erase sibling nodes")

    const auto lock = _mutex.makeScopedLock();
    std::stack<TIndex,std::vector<TIndex>> indexStack;
    indexStack.push(nodeIndex);

    while (!indexStack.empty()) {
        const auto i_node = indexStack.top();
        indexStack.pop();

        for (unsigned i_sibling=0; i_sibling<ContiguousTree::ChildrenPerNode; ++i_sibling) {
            Ref<TNode> r_node = this->_nodes[i_node + i_sibling];
            if (!r_node.isLeaf()) {
                indexStack.push(r_node.getChildBegin());
            }
            r_node.setNull();
        }

        _availableIndices.push(i_node);
    }

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor>
inline void
ContiguousTree<TNode,TIndex,TTag>::visit(const TFunctor& r_visitFunctor) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(!_nodes.empty())
    using IndexPair = std::pair<TIndex,TIndex>; // {i_node, depth}
    std::stack<IndexPair,DynamicArray<IndexPair>> indexStack;

    if (!_nodes.front().isLeaf()) {
        indexStack.emplace(0, 0);
    }

    while (!indexStack.empty()) {
        const auto indexPair = indexStack.top();
        indexStack.pop();

        for (TIndex i_sibling=0; i_sibling<TNode::ChildrenPerNode; ++i_sibling) {
            const TIndex i_node = indexPair.first + i_sibling;
            Ref<const Node> r_node = _nodes[i_node];
            if (!r_visitFunctor(r_node, indexPair.second)) {
                return;
            }
            if (!r_node.isLeaf()) {
                indexStack.emplace(r_node.getChildBegin(), indexPair.second + 1);
            }
        }
    } // while indexStack

    CIE_END_EXCEPTION_TRACING
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::FunctionWithSignature<bool, const TNode&, TIndex> TFunctor, concepts::ThreadPool TPool>
inline void
ContiguousTree<TNode,TIndex,TTag>::visit(const TFunctor& r_visitFunctor, TPool& r_pool) const
{
    CIE_THROW(NotImplementedException, "Multithreaded node visiting is not implemented yet")
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
template <concepts::Iterator<TNode> TIterator>
inline TIndex ContiguousTree<TNode,TIndex,TTag>::insert(TIterator it_nodeBegin)
{
    CIE_BEGIN_EXCEPTION_TRACING

    static_assert(TNode::ChildrenPerNode < std::numeric_limits<TIndex>::max(),
                  "Children cannot be indexed by the specified index type!");

    const auto lock = _mutex.makeScopedLock();

    TIndex i_begin;

    if (_availableIndices.empty()) {
        // Check for overflows
        const auto size = _nodes.size();
        i_begin = size;
        const auto extendedSize = _nodes.size() + TNode::ChildrenPerNode;
        CIE_CHECK(size < extendedSize && extendedSize <= std::numeric_limits<TIndex>::max(),
                  "Index overflow!")

        // _nodes.reserve(extendedSize); // <-- this isn't that relevant here
        for (TIndex i=0; i<ContiguousTree::ChildrenPerNode; ++i, ++it_nodeBegin) {
            _nodes.emplace_back(*it_nodeBegin);
        }
    } else {
        // Pop an available begin from the set.
        i_begin = _availableIndices.front();
        _availableIndices.pop();

        // Insert nodes
        auto it_node = _nodes.begin() + i_begin;
        const auto it_nodeEnd = it_node + ContiguousTree::ChildrenPerNode;
        for (; it_node!=it_nodeEnd; ++it_node, ++it_nodeBegin) {
            *it_node = *it_nodeBegin;
        }
    }

    return i_begin;

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
inline Size
ContiguousTree<TNode,TIndex,TTag>::size() const noexcept
{
    return std::count_if(_nodes.begin(),
                         _nodes.end(),
                         [](const auto& r_node){return !r_node.isNull();});
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline Bool
ContiguousTree<TNode,TIndex,TTag>::empty() const noexcept
{
    return this->size() == 0;
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline Size
ContiguousTree<TNode,TIndex,TTag>::numberOfNodes() const
{
    return _nodes.size();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline Size
ContiguousTree<TNode,TIndex,TTag>::numberOfNulls() const
{
    return _nodes.size() - this->size();
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline typename ContiguousTree<TNode,TIndex,TTag>::const_iterator
ContiguousTree<TNode,TIndex,TTag>::begin() const noexcept
{
    return Iterator(this->_nodes.begin(),
                    this->_nodes.begin(),
                    this->_nodes.end());
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline typename ContiguousTree<TNode,TIndex,TTag>::iterator
ContiguousTree<TNode,TIndex,TTag>::begin() noexcept
{
    return Iterator(this->_nodes.begin(),
                    this->_nodes.begin(),
                    this->_nodes.end());
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline typename ContiguousTree<TNode,TIndex,TTag>::const_iterator
ContiguousTree<TNode,TIndex,TTag>::end() const noexcept
{
    return Iterator(this->_nodes.end(),
                    this->_nodes.begin(),
                    this->_nodes.end());
}


template <concepts::ContiguousTreeNode TNode, concepts::Integer TIndex, class TTag>
inline typename ContiguousTree<TNode,TIndex,TTag>::iterator
ContiguousTree<TNode,TIndex,TTag>::end() noexcept
{
    return Iterator(this->_nodes.end(),
                    this->_nodes.begin(),
                    this->_nodes.end());
}


} // namespace cie::utils


#endif
