#ifdef CIE_ENABLE_SYCL

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/concurrency/inc/sycl.hpp"

// --- STL Includes ---
#include <iostream>


CIE_TEST_CASE("SYCL", "[concurrency]")
{
    {
        sycl::default_selector deviceSelector;
        sycl::queue queue(deviceSelector);

        std::cout << "platform\n"
                  << "\tname                : " << queue.get_device().get_platform().get_info<sycl::info::platform::name>() << "\n"
                  << "\tversion             : " << queue.get_device().get_platform().get_info<sycl::info::platform::version>() << "\n";

        std::cout << "device\n"
                  << "\tname                : " << queue.get_device().get_info<sycl::info::device::name>() << "\n"
                  << "\tmax work group size : " << queue.get_device().get_info<sycl::info::device::max_work_group_size>() << "\n"
                  << "\tmax compute units   : " << queue.get_device().get_info<sycl::info::device::max_compute_units>() << std::endl;
    }
}

#endif // CIE_ENABLE_SYCL
