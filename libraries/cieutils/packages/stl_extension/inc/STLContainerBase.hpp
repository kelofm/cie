#ifndef CIE_UTILS_STL_CONTAINER_BASE
#define CIE_UTILS_STL_CONTAINER_BASE

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/resize.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie {


///@addtogroup cieutils
///@{

template <concepts::Container TBase>
class STLContainerBase : protected TBase
{
public:
    using typename TBase::value_type;

    using typename TBase::size_type;

    using typename TBase::difference_type;

    using typename TBase::pointer;

    using typename TBase::const_pointer;

    using typename TBase::reference;

    using typename TBase::const_reference;

    using typename TBase::iterator;

    using typename TBase::const_iterator;

    using typename TBase::reverse_iterator;

    using typename TBase::const_reverse_iterator;

    CIE_DEFINE_CLASS_POINTERS(STLContainerBase)

public:
    ///@name Constructor / Destructor
    ///@{

    STLContainerBase(const std::initializer_list<value_type>& r_initializer)
    requires concepts::detail::HasPushBack<TBase, const value_type&, void>
    {
        CIE_BEGIN_EXCEPTION_TRACING
        utils::reserve(*this, r_initializer.size());
        std::copy(r_initializer.begin(), r_initializer.end(), std::back_inserter(*this));
        CIE_END_EXCEPTION_TRACING
    }

    STLContainerBase(const std::initializer_list<value_type>& r_initializer)
    requires (!concepts::detail::HasPushBack<TBase, const value_type&, void>)
    {
        CIE_BEGIN_EXCEPTION_TRACING
        utils::resize(*this, r_initializer.size());
        std::copy(r_initializer.begin(), r_initializer.end(), this->begin());
        CIE_END_EXCEPTION_TRACING
    }

    STLContainerBase(RightRef<TBase> r_rhs) noexcept
        : TBase(std::move(r_rhs))
    {}

    STLContainerBase(Ref<const TBase> r_rhs)
        : TBase(r_rhs)
    {}

    STLContainerBase(STLContainerBase<TBase>&& r_rhs) noexcept = default;

    STLContainerBase(const STLContainerBase<TBase>& r_rhs) = default;

    template <class ...TArgs>
    STLContainerBase(TArgs&&... r_arguments)
        : TBase(std::forward<TArgs>(r_arguments)...)
    {}

    STLContainerBase<TBase>& operator=(STLContainerBase<TBase>&& r_rhs) noexcept = default;

    STLContainerBase<TBase>& operator=(const STLContainerBase<TBase>& r_rhs) noexcept = default;

    virtual ~STLContainerBase()
    {}

    ///@}
    ///@name Common members
    ///@{

    size_type size() const noexcept
    {return TBase::size();}

    bool empty() const noexcept
    {return TBase::empty();}

    iterator begin() noexcept
    {return TBase::begin();}

    const_iterator begin() const noexcept
    {return TBase::begin();}

    const_iterator cbegin() const noexcept
    {return TBase::cbegin();}

    reverse_iterator rbegin() noexcept
    {return TBase::rbegin();}

    const_reverse_iterator rbegin() const noexcept
    {return TBase::rbegin();}

    const_reverse_iterator crbegin() const noexcept
    {return TBase::crbegin();}

    iterator end() noexcept
    {return TBase::end();}

    const_iterator end() const noexcept
    {return TBase::end();}

    const_iterator cend() const noexcept
    {return TBase::cend();}

    reverse_iterator rend() noexcept
    {return TBase::rend();}

    const_reverse_iterator rend() const noexcept
    {return TBase::rend();}

    const_reverse_iterator crend() const noexcept
    {return TBase::crend();}

    ///@}
    ///@name Optional members
    ///@{

    reference at(size_type index) noexcept
    requires concepts::detail::HasAt<TBase, size_type, reference>
    {return TBase::at(index);}

    const_reference at(size_type index) const noexcept
    requires concepts::detail::HasAt<const TBase, size_type, const_reference>
    {return TBase::at(index);}

    reference operator[](size_type index) noexcept
    requires concepts::detail::HasAccessOperator<TBase, size_type, reference>
    {return TBase::operator[](index);}

    const_reference operator[](size_type index) const noexcept
    requires concepts::detail::HasAccessOperator<const TBase, size_type, const_reference>
    {return TBase::operator[](index);}

    pointer data() noexcept
    requires concepts::detail::HasData<TBase>
    {return TBase::data();}

    const_pointer data() const noexcept
    requires concepts::detail::HasData<TBase>
    {return TBase::data();}

    reference front() noexcept
    requires concepts::detail::HasFront<TBase, reference>
    {return TBase::front();}

    const_reference front() const noexcept
    requires concepts::detail::HasFront<const TBase, const_reference>
    {return TBase::front();}

    reference back() noexcept
    requires concepts::detail::HasBack<TBase, reference>
    {return TBase::back();}

    const_reference back() const noexcept
    requires concepts::detail::HasBack<const TBase, const_reference>
    {return TBase::back();}

    void reserve(size_type capacity)
    requires concepts::detail::HasReserve<TBase>
    {TBase::reserve(capacity);}

    auto capacity() const noexcept
    requires concepts::detail::HasReserve<TBase>
    {return TBase::capacity();}

    void clear()
    requires concepts::detail::HasClear<TBase>
    {TBase::clear();}

    iterator erase(iterator position)
    requires concepts::detail::HasErase<TBase, iterator>
    {return TBase::erase(position);}

    iterator erase(iterator begin, iterator end)
    requires concepts::detail::HasErase<TBase, iterator>
    {return TBase::erase(begin, end);}

    void push_front(const value_type& r_value)
    requires concepts::detail::HasPushBack<TBase, const value_type&, void>
    {TBase::push_front(r_value);}

    void push_front(value_type&& r_value)
    requires concepts::detail::HasPushBack<TBase, const value_type&, void>
    {TBase::push_front(std::move(r_value));}

    void push_back(const value_type& r_value)
    requires concepts::detail::HasPushBack<TBase, const value_type&, void>
    {TBase::push_back(r_value);}

    void push_back(value_type&& r_value)
    requires concepts::detail::HasPushBack<TBase, const value_type&, void>
    {TBase::push_back(std::move(r_value));}

    template <class ...TArgs>
    reference emplace_front(TArgs&&... r_arguments)
    requires concepts::detail::HasEmplaceBack<TBase, TArgs...>
    {return TBase::emplace_back(std::forward<TArgs>(r_arguments)...);}

    template <class ...TArgs>
    reference emplace_back(TArgs&&... r_arguments)
    requires concepts::detail::HasEmplaceBack<TBase, TArgs...>
    {return TBase::emplace_back(std::forward<TArgs>(r_arguments)...);}

    void pop_front()
    requires concepts::detail::HasPopFront<TBase>
    {TBase::pop_front();}

    void pop_back()
    requires concepts::detail::HasPopBack<TBase>
    {TBase::pop_back();}

    void resize(size_type size)
    requires concepts::detail::HasResize<TBase, size_type>
    {TBase::resize(size);}

    //void swap(STLContainerBase<TBase>& r_rhs)
    //requires concepts::detail::HasSwap<TBase>
    //{TBase::swap(r_rhs);}

    ///@}

    friend bool operator<(const STLContainerBase& r_lhs, const STLContainerBase& r_rhs)
    {return static_cast<const TBase&>(r_lhs) < static_cast<const TBase&>(r_rhs);}

    friend bool operator!=(const STLContainerBase& r_lhs, const STLContainerBase& r_rhs)
    {return static_cast<const TBase&>(r_lhs) != static_cast<const TBase&>(r_rhs);}

private:
    explicit operator const TBase& () const {return *this;}

    explicit operator TBase& () {return *this;}
};

///@}

} // namespace cie


#endif