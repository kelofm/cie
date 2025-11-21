#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <array> // array
#include <iterator> // reverse_iterator


namespace cie {


template <class T, std::size_t TMaxSize>
class StackArray
{
public:
    using value_type                = T;
    using size_type                 = std::size_t;
    using difference_type           = std::ptrdiff_t;
    using pointer                   = Ptr<T>;
    using const_pointer             = Ptr<const T>;
    using reference                 = Ref<T>;
    using const_reference           = Ref<const T>;
    using iterator                  = pointer;
    using const_iterator            = const_pointer;
    using reverse_iterator          = std::reverse_iterator<iterator>;
    using const_reverse_iterator    = std::reverse_iterator<const_iterator>;

    constexpr StackArray() noexcept
    {}

    [[nodiscard]] constexpr size_type size() const noexcept {
        return _iEnd / sizeof(T);}

    [[nodiscard]] constexpr bool empty() const noexcept {
        return _iEnd == 0ul;}

    [[nodiscard]] constexpr size_type capacity() const noexcept {
        return TMaxSize;}

    [[nodiscard]] constexpr const_iterator begin() const noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data()));}

    [[nodiscard]] constexpr iterator begin() noexcept {
        return static_cast<iterator>(static_cast<Ptr<void>>(_data.data.data()));}

    [[nodiscard]] constexpr const_iterator cbegin() noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data()));}

    [[nodiscard]] constexpr const_iterator cbegin() const noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data()));}

    [[nodiscard]] constexpr const_reverse_iterator rbegin() const noexcept {
        return const_reverse_iterator(this->end());}

    [[nodiscard]] constexpr reverse_iterator rbegin() noexcept {
        return reverse_iterator(this->end());}

    [[nodiscard]] constexpr const_reverse_iterator crbegin() const noexcept {
        return const_reverse_iterator(this->end());}

    [[nodiscard]] constexpr const_reverse_iterator crbegin() noexcept {
        return const_reverse_iterator(this->end());}

    [[nodiscard]] constexpr const_iterator end() const noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data() + _iEnd));}

    [[nodiscard]] constexpr iterator end() noexcept {
        return static_cast<iterator>(static_cast<Ptr<void>>(_data.data.data() + _iEnd));}

    [[nodiscard]] constexpr const_iterator cend() noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data() + _iEnd));}

    [[nodiscard]] constexpr const_iterator cend() const noexcept {
        return static_cast<const_iterator>(static_cast<Ptr<const void>>(_data.data.data() + _iEnd));}

    [[nodiscard]] constexpr const_reverse_iterator rend() const noexcept {
        return const_reverse_iterator(this->begin());}

    [[nodiscard]] constexpr reverse_iterator rend() noexcept {
        return reverse_iterator(this->begin());}

    [[nodiscard]] constexpr const_reverse_iterator crend() const noexcept {
        return const_reverse_iterator(this->begin());}

    [[nodiscard]] constexpr const_reverse_iterator crend() noexcept {
        return const_reverse_iterator(this->begin());}

    [[nodiscard]] constexpr const_pointer data() const noexcept {
        return this->begin();}

    [[nodiscard]] constexpr pointer data() noexcept {
        return this->begin();}

    [[nodiscard]] constexpr const_reference front() const noexcept {
        return *this->begin();}

    [[nodiscard]] constexpr reference front() noexcept {
        return *this->begin();}

    [[nodiscard]] constexpr const_reference back() const noexcept {
        return *this->rbegin();}

    [[nodiscard]] constexpr reference back() noexcept {
        return *this->rbegin();}

    [[nodiscard]] constexpr const_reference at(size_type index) const noexcept {
        return *(this->begin() + index);}

    [[nodiscard]] constexpr reference at(size_type index) noexcept {
        return *(this->begin() + index);}

    [[nodiscard]] constexpr const_reference operator[](size_type index) const noexcept {
        return this->at(index);}

    [[nodiscard]] constexpr reference operator[](size_type index) noexcept {
        return this->at(index);}

    [[nodiscard]] constexpr bool reserve(size_type count) noexcept {
        if (count <= TMaxSize) {
            return true;
        } else {
            return false;
        }
    }

    [[nodiscard]] bool resize(size_type newSize) {
        const auto currentSize = this->size();
        if (currentSize < newSize) {
            if (this->reserve(newSize)) {
                // Call the default constructor of all the newly "allocated" items.
                if constexpr (!std::is_integral_v<T>) {
                    iterator itConstruct = this->begin() + currentSize;
                    const iterator itEnd = this->begin() + newSize + 1;
                    for (; itConstruct<itEnd; ++itConstruct) new (itConstruct) T();
                }
                _iEnd = newSize * sizeof(T);
            } else {
                return false;
            }
        } else if (newSize < currentSize) {
            // Call the destructor of erased objects.
            if constexpr (!std::is_integral_v<T>) {
                reverse_iterator itDestroy   = ++reverse_iterator(this->begin() + currentSize);
                const reverse_iterator itEnd = ++reverse_iterator(this->begin() + newSize + 1);
                for (; itDestroy < itEnd; ++itDestroy) itDestroy->~T();
            }
            _iEnd = newSize * sizeof(T);
        }
        return true;
    }

    template <class ...TArgs>
    [[nodiscard]] bool emplace_back(TArgs&&... rArgs) {
        const auto it = this->end();
        if (!this->reserve(this->size() + 1ul)) return false;
        new (it) T(std::forward<TArgs>(rArgs)...);
        _iEnd += sizeof(T);
    }

    [[nodiscard]] bool push_back(T&& r) {
        const auto it = this->end();
        if (!this->reserve(this->size() + 1ul)) return false;
        new (it) T(std::move(r));
        _iEnd += sizeof(T);
    }

    [[nodiscard]] bool push_back(const T& r) {
        const auto it = this->end();
        if (!this->reserve(this->size() + 1ul)) return false;
        new (it) T(r);
        _iEnd += sizeof(T);
    }

    void pop_back() {
        if constexpr (!std::is_integral_v<T>) this->back().~T();
        _iEnd -= sizeof(T);
    }

    void clear() {
        (void)this->resize(0ul);
    }

private:
    using ObjectArray = std::array<T,TMaxSize>;
    using ByteArray = std::array<std::byte,TMaxSize * sizeof(T)>;
    struct alignas(ByteArray) Data {
        ByteArray data;
    }; // struct Data

    Data _data;

    std::size_t _iEnd;
}; // class StackArray


} // namespace cie
