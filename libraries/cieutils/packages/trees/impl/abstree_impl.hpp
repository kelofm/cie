#ifndef CIE_CIEUTILS_ABS_TREE_IMPL_HPP
#define CIE_CIEUTILS_ABS_TREE_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/getReference.hpp"


namespace cie::utils {


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
AbsTree<TSelf,TContainer,TStored,TArgs...>::AbsTree(Size level) noexcept
    : _level(level)
{
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
AbsTree<TSelf,TContainer,TStored,TArgs...>::AbsTree() noexcept
    : _level(SIZE_MAX)
{
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline const typename AbsTree<TSelf,TContainer,TStored,TArgs...>::StoredContainer&
AbsTree<TSelf,TContainer,TStored,TArgs...>::children() const
{
    return _children;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline typename AbsTree<TSelf,TContainer,TStored,TArgs...>::StoredContainer&
AbsTree<TSelf,TContainer,TStored,TArgs...>::children()
{
    return this->_children;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline const TSelf&
AbsTree<TSelf,TContainer,TStored,TArgs...>::child(Size index) const
{
    const auto& reference = getRef(_children[index]);
    return reference;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline TSelf&
AbsTree<TSelf,TContainer,TStored,TArgs...>::child(Size index)
{
    auto& reference = getRef(_children[index]);
    return reference;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline Size
AbsTree<TSelf,TContainer,TStored,TArgs...>::level() const
{
    return _level;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor>
inline bool
AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& r_visitor)
{
    bool result = r_visitor(dynamic_cast<TSelf*>(this));

    if (result)
        for (auto& r_child : this->_children)
        {
            result = getRef(r_child).visit(r_visitor);
            if (!result)
                break;
        }

    return result;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor>
inline bool
AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& r_visitor) const
{
    bool result = r_visitor(dynamic_cast<const TSelf*>(this));

    if (result)
        for (auto& r_child : this->_children)
        {
            result = getRef(r_child).visit(r_visitor);
            if (!result)
                break;
        }

    return result;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor, class TPool>
inline void
AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& r_visitor, TPool& r_threadPool)
{
    r_threadPool.queueJob(std::bind(r_visitor, this));
    for (auto& r_child : _children)
        getRef(r_child).visit(r_visitor);
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor, class TPool>
inline void
AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& r_visitor, TPool& r_threadPool) const
{
    r_threadPool.queueJob(std::bind(r_visitor, this));
    for (const auto& r_child : _children)
        getRef(r_child).visit(r_visitor);
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline void
AbsTree<TSelf,TContainer,TStored,TArgs...>::clear()
{
    _children.clear();
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class ...Ts>
requires concepts::Pointer<TStored>
inline TSelf&
AbsTree<TSelf,TContainer,TStored,TArgs...>::emplaceChild(Ts&&... r_arguments)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _children.emplace_back(new TSelf(std::forward<Ts>(r_arguments)...));
    return *_children.back();

    CIE_END_EXCEPTION_TRACING
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class ...Ts>
requires concepts::NonPointer<TStored>
inline TSelf&
AbsTree<TSelf,TContainer,TStored,TArgs...>::emplaceChild(Ts&&... r_arguments)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _children.emplace_back(std::forward<Ts>(r_arguments)...);
    return *_children.back();

    CIE_END_EXCEPTION_TRACING
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
inline bool
AbsTree<TSelf,TContainer,TStored,TArgs...>::isLeaf() const noexcept
{
    if constexpr (concepts::Pointer<TStored>) {
        // Has no children, or all children are nullptrs

        if (this->_children.empty())
            return true;

        for (const auto& r_child : this->_children)
            if (r_child)
                return false;

        return true;
    } else {
        return _children.empty();
    }
}


} // namespace cie::utils

#endif