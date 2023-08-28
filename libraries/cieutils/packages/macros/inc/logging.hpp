#ifndef CIE_UTILS_MACROS_LOGGING_HPP
#define CIE_UTILS_MACROS_LOGGING_HPP

// --- Utility Includes ---
#ifdef CIE_ENABLE_PROFILING
#include "packages/logging/inc/Profiler.hpp"
#endif


#ifdef CIE_ENABLE_PROFILING
    #define CIE_DEFINE_SCOPE_PROFILER_MACRO(CIE_TIME_UNIT)                                                                  \
        thread_local static auto& CIE_SCOPE_PROFILED_ITEM =                                                                 \
            cie::utils::ProfilerSingleton<CIE_TIME_UNIT>::get().create(std::source_location::current());                    \
        auto CIE_SCOPE_PROFILER = cie::utils::ProfilerSingleton<CIE_TIME_UNIT>::get().profile(CIE_SCOPE_PROFILED_ITEM);

    #define CIE_PROFILE_SCOPE_MILLI CIE_DEFINE_SCOPE_PROFILER_MACRO(std::chrono::milliseconds)

    #define CIE_PROFILE_SCOPE_MICRO CIE_DEFINE_SCOPE_PROFILER_MACRO(std::chrono::microseconds)

    #define CIE_PROFILE_SCOPE_NANO CIE_DEFINE_SCOPE_PROFILER_MACRO(std::chrono::nanoseconds)

    #define CIE_PROFILE_SCOPE CIE_PROFILE_SCOPE_MICRO

    //#undef CIE_DEFINE_SCOPE_PROFILER_MACRO

#else
    #define CIE_PROFILE_SCOPE_MILLI

    #define CIE_PROFILE_SCOPE_MICRO

    #define CIE_PROFILE_SCOPE_NANO

    #define CIE_PROFILE_SCOPE
#endif


#endif
