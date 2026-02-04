#pragma once

// --- External Includes ---
#ifdef CIE_ENABLE_SYCL
    #ifndef SYCL_LANGUAGE_VERSION
        #define SYCL_LANGUAGE_VERSION 2020
    #endif
    #include <sycl/sycl.hpp>
#endif

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <memory>
#include <format>


namespace cie {


template <class T>
using DynamicSharedArray =
#ifdef CIE_ENABLE_SYCL
    DynamicArray<T,sycl::usm_allocator<T,sycl::usm::alloc::shared>>
#else
    DynamicArray<T>
#endif
;


#ifdef CIE_ENABLE_SYCL
template <class T>
class DeviceMemoryDeleter {
public:
    constexpr DeviceMemoryDeleter(sycl::queue& rQueue) noexcept
        : _pQueue(&rQueue)
    {}

    constexpr DeviceMemoryDeleter() noexcept
        : _pQueue(nullptr)
    {}

    void operator()(T* p) const {
        //utils::LoggerSingleton::get().log(std::format(
        //    "SYCL: deallocate {}",
        //    static_cast<void*>(p)));
        sycl::free(p, *_pQueue);
    }

private:
    Ptr<sycl::queue> _pQueue;
}; // class DeviceMemoryDeleter


template <class T>
using DeviceMemory = std::unique_ptr<
    T,
    DeviceMemoryDeleter<T>>;

template <class T>
DeviceMemory<T> makeDeviceMemory(std::size_t size, Ref<sycl::queue> rQueue) {
    T* p = sycl::malloc_device<T>(size, rQueue);
    //utils::LoggerSingleton::get().log(std::format(
    //    "SYCL: allocate {} bytes at {}",
    //    size * sizeof(T),
    //    static_cast<void*>(p)));
    return DeviceMemory<T>(p, DeviceMemoryDeleter<T>(rQueue));
}
#endif


} // namespace cie
