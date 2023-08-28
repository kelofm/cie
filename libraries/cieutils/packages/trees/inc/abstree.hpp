#ifndef CIE_CIEUTILS_ABS_TREE_HPP
#define CIE_CIEUTILS_ABS_TREE_HPP

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


// Tree base class
template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
class AbsTree
{
public:
    using Self = TSelf;

    using Stored = TStored;

    using StoredContainer = TContainer<TStored,TArgs...>;

    //CIE_DEFINE_CLASS_POINTERS(TSelf)

public:

    AbsTree();

    AbsTree(Size level);

    AbsTree(AbsTree&& r_rhs) = default;

    AbsTree(const AbsTree& r_rhs) = default;

    AbsTree& operator=(AbsTree&& r_rhs) = default;

    AbsTree& operator=(const AbsTree& r_rhs) = default;

    virtual ~AbsTree() = default;

    /** Send a function down the tree and execute it on all nodes while it returns true
     *  @param r_visitor lambda taking a node pointer and returning a bool
     */
    template <class TVisitor>
    bool visit(TVisitor&& r_visitor);

    /** Send a function down the tree and execute it on all nodes while it returns true
     *  @param r_visitor lambda taking a const node pointer and returning a bool
     */
    template <class TVisitor>
    bool visit(TVisitor&& r_visitor) const;

    /** Send a function down the tree and execute it on all nodes
     *  @param r_visitor lambda taking a const node pointer and returning a bool
     *  @param r_threadPool thread pool to assign the evaluations to
     */
    template <class TVisitor, class TPool>
    void visit(TVisitor&& r_visitor, TPool& r_threadPool);

    /** Send a function down the tree and execute it on all nodes
     *  @param r_visitor lambda taking a const node pointer and returning a bool
     *  @param r_threadPool thread pool to assign the evaluations to
     */
    template <class TVisitor, class TPool>
    void visit(TVisitor&& r_visitor, TPool& r_threadPool) const;

    /// Clear children (non-recursive)
    virtual void clear();

    /// Check whether this node is a leaf node
    bool isLeaf() const
    requires concepts::Pointer<TStored>;

    /// Check whether this node is a leaf node
    bool isLeaf() const
    requires concepts::NonPointer<TStored>;

    template <class ...Ts>
    requires concepts::Pointer<TStored>
    TSelf& emplaceChild(Ts&&... r_arguments);

    template <class ...Ts>
    requires concepts::NonPointer<TStored>
    TSelf& emplaceChild(Ts&&... r_arguments);

    Size level() const;

    const StoredContainer& children() const;

    StoredContainer& children();

    const TSelf& child(Size index) const;

    TSelf& child(Size index);

protected:
    StoredContainer _children;

    Size _level;
};


} // namespace cie::utils

#include "packages/trees/impl/abstree_impl.hpp"

#endif