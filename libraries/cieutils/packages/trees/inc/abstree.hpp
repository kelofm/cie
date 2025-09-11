#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <functional>
#include <memory>

namespace cie::utils {


/// Basic tree class.
/// @ingroup cieutils
template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
class AbsTree
{
public:
    using Self = TSelf;

    using Stored = TStored;

    using StoredContainer = TContainer<TStored,TArgs...>;

    //CIE_DEFINE_CLASS_POINTERS(TSelf)

public:
    AbsTree() noexcept;

    AbsTree(Size level) noexcept;

    AbsTree(AbsTree&& r_rhs) noexcept = default;

    AbsTree(const AbsTree& r_rhs) = default;

    AbsTree& operator=(AbsTree&& r_rhs) noexcept = default;

    AbsTree& operator=(const AbsTree& r_rhs) = default;

    virtual ~AbsTree() = default;

    /** Send a function down the tree and execute it on all nodes while it returns true
     *  @param rVisitor lambda taking a node pointer and returning a bool
     */
    template <class TVisitor>
    bool visit(TVisitor&& rVisitor);

    /** Send a function down the tree and execute it on all nodes while it returns true
     *  @param rVisitor lambda taking a const node pointer and returning a bool
     */
    template <class TVisitor>
    bool visit(TVisitor&& rVisitor) const;

    /** Send a function down the tree and execute it on all nodes
     *  @param rVisitor lambda taking a const node pointer and returning a bool
     *  @param rThreadPool thread pool to assign the evaluations to
     */
    template <class TVisitor, class TPool>
    void visit(TVisitor&& rVisitor, TPool& rThreadPool);

    /** Send a function down the tree and execute it on all nodes
     *  @param rVisitor lambda taking a const node pointer and returning a bool
     *  @param rThreadPool thread pool to assign the evaluations to
     */
    template <class TVisitor, class TPool>
    void visit(TVisitor&& rVisitor, TPool& rThreadPool) const;

    /// Clear children (non-recursive)
    virtual void clear();

    /// Check whether this node is a leaf node.
    bool isLeaf() const noexcept;

    template <class ...Ts>
    requires concepts::Pointer<TStored>
    TSelf& emplaceChild(Ts&&... rArguments);

    template <class ...Ts>
    requires concepts::NonPointer<TStored>
    TSelf& emplaceChild(Ts&&... rArguments);

    Size level() const;

    const StoredContainer& children() const;

    StoredContainer& children();

    const TSelf& child(Size index) const;

    TSelf& child(Size index);

protected:
    StoredContainer _children;

    std::uint16_t _level;
};


} // namespace cie::utils

#include "packages/trees/impl/abstree_impl.hpp"
