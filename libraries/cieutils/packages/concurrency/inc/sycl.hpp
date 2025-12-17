#pragma once

#ifdef CIE_ENABLE_SYCL
    #ifndef SYCL_LANGUAGE_VERSION
        #define SYCL_LANGUAGE_VERSION 2020
    #endif
    #include <sycl/sycl.hpp>
#endif
