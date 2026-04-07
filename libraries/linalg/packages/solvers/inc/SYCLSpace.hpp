#pragma once
#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/concurrency/inc/sycl.hpp"


namespace cie::linalg {



template <class T>
class SYCLSpace;



template <class T>
class SYCLVector {
public:
    SYCLVector() noexcept = default;

    SYCLVector(
        DeviceMemory<T>&& rMemory,
        std::size_t size);

    [[nodiscard]] std::size_t size() const noexcept;

    [[nodiscard]] Ptr<const T> get() const noexcept;

    [[nodiscard]] Ptr<T> get() noexcept;

private:
    DeviceMemory<T> _pMemory;

    std::size_t _size;
}; // class SYCLVector



template <class T>
class SYCLView {
public:
    SYCLView() noexcept = default;

    SYCLView(
        std::conditional_t<
            std::is_const_v<T>,
            Ref<const SYCLVector<std::remove_const_t<T>>>,
            Ref<SYCLVector<std::remove_const_t<T>>>
        > rVector) noexcept;

    [[nodiscard]] std::size_t size() const noexcept;

    [[nodiscard]] Ptr<const T> get() const noexcept;

    [[nodiscard]] Ptr<T> get() noexcept
    requires (!std::is_const_v<T>);

    [[nodiscard]] operator SYCLView<const T> () const noexcept;

private:
    std::conditional_t<
        std::is_const_v<T>,
        Ptr<const SYCLVector<std::remove_const_t<T>>>,
        Ptr<SYCLVector<std::remove_const_t<T>>>
    > _pVector;
}; // class SYCLView



template <class T>
class SYCLSpace {
public:
    using Vector = SYCLVector<T>;

    using Value = T;

    using VectorView = SYCLView<T>;

    using ConstVectorView = SYCLView<const T>;

    SYCLSpace() noexcept = default;

    SYCLSpace(std::shared_ptr<sycl::queue> pQueue);

    [[nodiscard]] static VectorView view(Ref<Vector> rVector) noexcept;

    [[nodiscard]] static ConstVectorView view(Ref<const Vector> rVector) noexcept;

    [[nodiscard]] static std::size_t size(ConstVectorView view) noexcept;

    [[nodiscard]] Vector makeVector(std::size_t size) const;

    [[nodiscard]] Value innerProduct(
        ConstVectorView left,
        ConstVectorView right) const;

    void scale(
        VectorView target,
        ConstVectorView source,
        Value scale) const;

    void scale(
        VectorView view,
        Value value) const;

    void add(
        VectorView target,
        ConstVectorView source,
        Value scale) const;

    void assign(
        VectorView target,
        ConstVectorView source) const;

    void assign(
        VectorView target,
        std::span<const T> source) const;

    void assign(
        std::span<T> target,
        ConstVectorView source) const;

    void fill(
        VectorView target,
        Value value) const;

    std::shared_ptr<sycl::queue> getQueue() const;

private:
    std::shared_ptr<sycl::queue> _pQueue;

    std::size_t _maxWorkGroupSize;
}; // class SYCLSpace


} // namespace cie::linalg

#endif
