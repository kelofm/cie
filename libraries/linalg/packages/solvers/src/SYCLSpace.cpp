#ifdef CIE_ENABLE_SYCL

// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLSpace.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <format>


namespace cie::linalg {


template <class T>
SYCLVector<T>::SYCLVector(
    DeviceMemory<T>&& rMemory,
    std::size_t size)
        :   _pMemory(std::move(rMemory)),
            _size(size)
{}


template <class T>
std::size_t SYCLVector<T>::size() const noexcept {
    return _size;
}


template <class T>
Ptr<const T> SYCLVector<T>::get() const noexcept {
    return _pMemory.get();
}


template <class T>
Ptr<T> SYCLVector<T>::get() noexcept {
    return _pMemory.get();
}


template <class T>
SYCLView<T>::SYCLView(
    std::conditional_t<
        std::is_const_v<T>,
        Ref<const SYCLVector<std::remove_const_t<T>>>,
        Ref<SYCLVector<std::remove_const_t<T>>>
    > rVector) noexcept
        : _pVector(&rVector)
{}


template <class T>
std::size_t SYCLView<T>::size() const noexcept {
    return _pVector->size();
}


template <class T>
Ptr<const T> SYCLView<T>::get() const noexcept {
    return _pVector->get();
}


template <class T>
Ptr<T> SYCLView<T>::get() noexcept
requires (!std::is_const_v<T>) {
    return _pVector->get();
}


template <class T>
SYCLView<T>::operator SYCLView<const T> () const noexcept {
    return SYCLView<const T>(*_pVector);
}


template <class T>
SYCLSpace<T>::SYCLSpace(std::shared_ptr<sycl::queue> pQueue)
    :   _pQueue(pQueue),
        _maxWorkGroupSize(pQueue->get_device().get_info<sycl::info::device::max_work_group_size>())
{}


template <class T>
typename SYCLSpace<T>::VectorView SYCLSpace<T>::view(Ref<Vector> rVector) noexcept {
    return VectorView(rVector);
}


template <class T>
typename SYCLSpace<T>::ConstVectorView SYCLSpace<T>::view(Ref<const Vector> rVector) noexcept {
    return ConstVectorView(rVector);
}


template <class T>
std::size_t SYCLSpace<T>::size(ConstVectorView view) noexcept {
    return view.size();
}


template <class T>
typename SYCLSpace<T>::Vector SYCLSpace<T>::makeVector(std::size_t size) const {
    return Vector(
        makeDeviceMemory<T>(
            size,
            *_pQueue),
        size);
}


template <class T>
typename SYCLSpace<T>::Value SYCLSpace<T>::innerProduct(
    ConstVectorView left,
    ConstVectorView right) const {
        CIE_CHECK(
            left.size() == right.size(),
            std::format(
                "incompatible array sizes for dot product: {} != {}",
                left.size(), right.size()))

        const std::size_t componentCount = left.size();
        if (componentCount == 0) return 0;

        const std::size_t componentsPerItem = std::min<std::size_t>(
            0x10,
            componentCount);
        const std::size_t workItemCount = (componentCount + componentsPerItem - 1) / componentsPerItem;

        CIE_BEGIN_EXCEPTION_TRACING
            const sycl::range<1> range(workItemCount);
            Ptr<const T> pLeftBegin = left.get();
            Ptr<const T> pRightBegin = right.get();

            // Initialize the result on the device.
            auto deviceResult = this->makeVector(1);
            this->fill(deviceResult, 0);

            // Compute the inner product.
            _pQueue->submit([&] (Ref<sycl::handler> rHandler) {
                auto reduction = sycl::reduction(
                    deviceResult.get(),
                    sycl::plus<T>());
                rHandler.parallel_for(
                    range,
                    reduction,
                    [componentsPerItem, componentCount, pLeftBegin, pRightBegin] (
                        sycl::item<1> it,
                        auto& rReduction) {
                            const std::size_t iItem = it.get_linear_id();
                            const std::size_t iComponentBegin = std::min<std::size_t>(
                                iItem * componentsPerItem,
                                componentCount);
                            const std::size_t iComponentEnd = std::min<std::size_t>(
                                iComponentBegin + componentsPerItem,
                                componentCount);
                            const T contribution = std::inner_product(
                                pLeftBegin + iComponentBegin,
                                pLeftBegin + iComponentEnd,
                                pRightBegin + iComponentBegin,
                                static_cast<T>(0));
                            rReduction += contribution;
                        });
            }).wait_and_throw();

            // Fetch the result from the device.
            T result = 0;
            this->assign({&result, 1}, deviceResult);
            return result;
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::scale(
    VectorView target,
    ConstVectorView source,
    Value scale) const {
        CIE_CHECK(target.size() == source.size(), "")
        if (!target.size()) return;

        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<T> pTargetBegin = target.get();
            Ptr<const T> pSourceBegin = source.get();
            const std::size_t size = target.size();
            auto job = [pTargetBegin, pSourceBegin, size, this] (const auto& op) -> void {
                _pQueue->parallel_for(
                    sycl::range<1>(size),
                    [pTargetBegin, pSourceBegin, op] (sycl::item<1> it) {
                        const std::size_t iComponent = it.get_linear_id();
                        pTargetBegin[iComponent] = op(pTargetBegin[iComponent], pSourceBegin[iComponent]);
                    }).wait_and_throw();
            }; // job

            if (scale == static_cast<T>(1)) {
                job([] (T t, T s) -> T {return t * s;});
            } else if (scale == static_cast<T>(-1)) {
                job([] (T t, T s) {return -t * s;});
            } else if (scale == static_cast<T>(0)) {
                job([] (T, T) {return 0;});
            } else {
                job([scale] (T t, T s) {return scale * t * s;});
            }
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::scale(
    VectorView target,
    Value scale) const {
        if (!target.size()) return;

        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<T> pTargetBegin = target.get();
            const std::size_t size = target.size();
            auto job = [pTargetBegin, size, this] (const auto& op) -> void {
                _pQueue->parallel_for(
                    sycl::range<1>(size),
                    [pTargetBegin, op] (sycl::item<1> it) {
                        const std::size_t iComponent = it.get_linear_id();
                        pTargetBegin[iComponent] = op(pTargetBegin[iComponent]);
                    }).wait_and_throw();
            }; // job

            if (scale == static_cast<T>(1)) {
                return;
            } else if (scale == static_cast<T>(-1)) {
                job([] (T target) -> T {return -target;});
            } else if (scale == static_cast<T>(0)) {
                job([] (T) -> T {return 0;});
            } else {
                job([scale] (T target) -> T {return target * scale;});
            }
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::add(
    VectorView target,
    ConstVectorView source,
    Value scale) const {
        CIE_CHECK(target.size() == source.size(), "")
        if (!target.size()) return;

        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<T> pTargetBegin = target.get();
            Ptr<const T> pSourceBegin = source.get();
            const std::size_t size = target.size();
            auto job = [pTargetBegin, pSourceBegin, size, this] (const auto& op) -> void {
                _pQueue->parallel_for(
                    sycl::range<1>(size),
                    [pTargetBegin, pSourceBegin, op] (sycl::item<1> it) {
                        const std::size_t iComponent = it.get_linear_id();
                        op(pTargetBegin[iComponent], pSourceBegin[iComponent]);
                    }).wait_and_throw();
            }; // job

            if (scale == static_cast<T>(1)) {
                job([] (Ref<T> rTarget, T source) {rTarget += source;});
            } else if (scale == static_cast<T>(-1)) {
                job([] (Ref<T> rTarget, T source) {rTarget -= source;});
            } else if (scale == static_cast<T>(0)) {
                job([] (Ref<T> rTarget, T source) {rTarget = source;});
            } else {
                job([scale] (Ref<T> rTarget, T source) {rTarget += scale * source;});
            }
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::fill(
    VectorView target,
    Value value) const {
        if (!target.size()) return;
        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<T> pTargetBegin = target.get();
            const std::size_t size = target.size();
            _pQueue->parallel_for(
                sycl::range<1>(size),
                [pTargetBegin, value] (sycl::item<1> it) -> void {
                    pTargetBegin[it.get_linear_id()] = value;
                }).wait_and_throw();
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::assign(
    VectorView target,
    ConstVectorView source) const {
        CIE_CHECK(target.size() == source.size(), "")
        if (!target.size()) return;
        CIE_BEGIN_EXCEPTION_TRACING
            Ptr<T> pTargetBegin = target.get();
            Ptr<const T> pSourceBegin = source.get();
            const std::size_t size = target.size();
            if (pTargetBegin == pSourceBegin || !size) return;

            _pQueue->submit([pTargetBegin, pSourceBegin, size] (Ref<sycl::handler> rHandler) {
                rHandler.parallel_for(
                    sycl::range<1>(size),
                    [pTargetBegin, pSourceBegin] (sycl::item<1> it) {
                        const std::size_t iComponent = it.get_linear_id();
                        pTargetBegin[iComponent] = pSourceBegin[iComponent];
                    });
            }).wait_and_throw();
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::assign(
    VectorView target,
    std::span<const T> source) const {
        CIE_CHECK(target.size() == source.size(), "")
        CIE_BEGIN_EXCEPTION_TRACING
            _pQueue->copy(
                source.data(),
                target.get(),
                source.size()
            ).wait_and_throw();
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void SYCLSpace<T>::assign(
    std::span<T> target,
    ConstVectorView source) const {
        CIE_CHECK(target.size() == source.size(), "")
        CIE_BEGIN_EXCEPTION_TRACING
            _pQueue->copy(
                source.get(),
                target.data(),
                source.size()
            ).wait_and_throw();
        CIE_END_EXCEPTION_TRACING
}


template <class T>
std::shared_ptr<sycl::queue> SYCLSpace<T>::getQueue() const {
    return _pQueue;
}


#define CIE_DEFINE_SYCL_SPACE(T)        \
    template class SYCLSpace<T>;        \
    template class SYCLVector<T>;       \
    template class SYCLView<T>;         \
    template class SYCLView<const T>;


CIE_DEFINE_SYCL_SPACE(std::uint16_t)

CIE_DEFINE_SYCL_SPACE(int)

CIE_DEFINE_SYCL_SPACE(std::size_t)

CIE_DEFINE_SYCL_SPACE(float)

CIE_DEFINE_SYCL_SPACE(double)


#undef CIE_DEFINE_SYCL_SPACE


} // namespace cie::linalg

#endif
