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


namespace cie {


template <class T>
using DynamicSharedArray =
#ifdef CIE_ENABLE_SYCL
    DynamicArray<T,sycl::usm_allocator<T,sycl::usm::alloc::shared>>
#else
    DynamicArray<T>
#endif
;


} // namespace cie
