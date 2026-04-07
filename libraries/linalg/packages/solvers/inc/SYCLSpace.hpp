#pragma once
#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/concurrency/inc/sycl.hpp"


namespace cie::linalg {


template <class T>
class SYCLSpace {
public:
    class Vector;

private:
    template <class TT>
    class ViewBase {
    public:
        ViewBase() noexcept = default;

        ViewBase(
            std::conditional_t<
                std::is_const_v<TT>,
                Ref<const Vector>,
                Ref<Vector>
            > rVector) noexcept;

        [[nodiscard]] std::size_t size() const noexcept;

        [[nodiscard]] Ptr<const T> get() const noexcept;

        [[nodiscard]] Ptr<T> get() noexcept
        requires (!std::is_const_v<TT>);

    private:
        std::conditional_t<
            std::is_const_v<TT>,
            Ptr<const Vector>,
            Ptr<Vector>
        > _pVector;
    }; // class ViewBase
public:
    class Vector {
    public:
        Vector() noexcept = default;

    private:
        friend class SYCLSpace<T>;

        Vector(
            DeviceMemory<T>&& rMemory,
            std::size_t size);

        [[nodiscard]] std::size_t size() const noexcept;

        [[nodiscard]] Ptr<const T> get() const noexcept;

        [[nodiscard]] Ptr<T> get() noexcept;

        DeviceMemory<T> _pMemory;

        std::size_t _size;
    }; // class Vector

    using Value = T;

    using VectorView = ViewBase<T>;

    using ConstVectorView = ViewBase<const T>;

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
