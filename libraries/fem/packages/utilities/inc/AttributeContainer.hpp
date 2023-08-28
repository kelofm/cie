#ifndef CIE_FEM_ATTRIBUTE_CONTAINER_HPP
#define CIE_FEM_ATTRIBUTE_CONTAINER_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <tuple>
#include <utility>


namespace cie::fem {


template <class ...TValues>
class AttributeContainer
{
private:
    using Tuple = std::tuple<DynamicArray<TValues>...>;

    constexpr static inline auto TupleIndices = std::make_index_sequence<std::tuple_size_v<Tuple>>();

public:
    using Values = std::tuple<TValues...>;

public:
    AttributeContainer() noexcept = default;

    AttributeContainer(Size size);

    void resize(Size size);

    void erase(Size i_begin, Size i_end);

    void erase(Size index);

    void clear();

    void reserve(Size capacity);

    void push_back(TValues... values);

    Size size() const noexcept;

    Size capacity() const noexcept;

    bool empty() const noexcept;

    std::tuple<TValues...> get(Size index) const noexcept;

    template <class TValue>
    TValue at(Size index) const noexcept;

    template <class TValue>
    Ref<TValue> at(Size index) noexcept;

    template <class ...TArguments, class TVisitor>
    auto visit(Size index, TVisitor&& r_visitor) const
    {
        return r_visitor(std::get<DynamicArray<TArguments>>(_containers)[index]...);
    }

    template <class TValue>
    typename DynamicArray<TValue>::const_iterator begin() const noexcept;

    template <class TValue>
    typename DynamicArray<TValue>::iterator begin() noexcept;

    template <class TValue>
    typename DynamicArray<TValue>::const_iterator end() const noexcept;

    template <class TValue>
    typename DynamicArray<TValue>::iterator end() noexcept;

private:
    Tuple _containers;
}; // class AttributeContainer


} // namespace cie::fem




namespace cie::concepts {


namespace impl {
template <class T>
struct AttributeContainer
{static constexpr bool value = false;};

template <class ...TValues>
struct AttributeContainer<fem::AttributeContainer<TValues...>>
{static constexpr bool value = true;};
} // namespace impl


template <class T>
concept AttributeContainer = impl::AttributeContainer<T>::value;


} // namespace cie::concepts



#include "packages/utilities/impl/AttributeContainer_impl.hpp"

#endif
