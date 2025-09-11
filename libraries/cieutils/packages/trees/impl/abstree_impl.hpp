#pragma once

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
const typename AbsTree<TSelf,TContainer,TStored,TArgs...>::StoredContainer&
AbsTree<TSelf,TContainer,TStored,TArgs...>::children() const
{
    return _children;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
typename AbsTree<TSelf,TContainer,TStored,TArgs...>::StoredContainer&
AbsTree<TSelf,TContainer,TStored,TArgs...>::children()
{
    return this->_children;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
const TSelf& AbsTree<TSelf,TContainer,TStored,TArgs...>::child(Size index) const
{
    const auto& reference = getRef(_children[index]);
    return reference;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
TSelf& AbsTree<TSelf,TContainer,TStored,TArgs...>::child(Size index)
{
    auto& reference = getRef(_children[index]);
    return reference;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
Size AbsTree<TSelf,TContainer,TStored,TArgs...>::level() const
{
    return _level;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor>
bool AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& rVisitor)
{
    bool result = rVisitor(dynamic_cast<TSelf*>(this));

    if (result)
        for (auto& rChild : this->_children) {
            getRef(rChild).visit(rVisitor);
        }

    return result;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor>
bool AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& rVisitor) const
{
    bool result = rVisitor(dynamic_cast<const TSelf*>(this));

    if (result)
        for (auto& rChild : this->_children) {
            getRef(rChild).visit(rVisitor);
        }

    return result;
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor, class TPool>
void AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& rVisitor, TPool& rThreadPool)
{
    rThreadPool.queueJob(std::bind(rVisitor, this));
    for (auto& rChild : _children)
        getRef(rChild).visit(rVisitor);
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class TVisitor, class TPool>
void AbsTree<TSelf,TContainer,TStored,TArgs...>::visit(TVisitor&& rVisitor, TPool& rThreadPool) const
{
    rThreadPool.queueJob(std::bind(rVisitor, this));
    for (const auto& rChild : _children)
        getRef(rChild).visit(rVisitor);
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
void AbsTree<TSelf,TContainer,TStored,TArgs...>::clear()
{
    _children.clear();
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class ...Ts>
requires concepts::Pointer<TStored>
TSelf& AbsTree<TSelf,TContainer,TStored,TArgs...>::emplaceChild(Ts&&... rArguments)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _children.emplace_back(new TSelf(std::forward<Ts>(rArguments)...));
    return *_children.back();

    CIE_END_EXCEPTION_TRACING
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
template <class ...Ts>
requires concepts::NonPointer<TStored>
TSelf& AbsTree<TSelf,TContainer,TStored,TArgs...>::emplaceChild(Ts&&... rArguments)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _children.emplace_back(std::forward<Ts>(rArguments)...);
    return *_children.back();

    CIE_END_EXCEPTION_TRACING
}


template <class TSelf, template <class ...> class TContainer, class TStored, class ...TArgs>
bool AbsTree<TSelf,TContainer,TStored,TArgs...>::isLeaf() const noexcept
{
    if constexpr (concepts::Pointer<TStored>) {
        // Has no children, or all children are nullptrs

        if (this->_children.empty())
            return true;

        for (const auto& rChild : this->_children)
            if (rChild)
                return false;

        return true;
    } else {
        return _children.empty();
    }
}


} // namespace cie::utils
