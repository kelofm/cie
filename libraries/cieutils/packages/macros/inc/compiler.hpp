#ifndef CIE_UTILS_COMPILER_MACROS_HPP
#define CIE_UTILS_COMPILER_MACROS_HPP


#define CIE_STRINGIFY(ARG) #ARG


#if defined(__GNUC__)
    #define CIE_BEGIN_DISABLE_COMPILER_WARNING(NAME) _Pragma(CIE_STRINGIFY(GCC diagnostic push ignored NAME))
    #define CIE_END_DISABLE_COMPILER_WARNING _Pragma("GCC diagnostic pop")
#elif defined(__clang__)
    #define CIE_BEGIN_DISABLE_COMPILER_WARNING(NAME) _Pragma(CIE_STRINGIFY(clang diagnostic push NAME))
    #define CIE_END_DISABLE_COMPILER_WARNING _Pragma("clang diagnostic pop")
#endif


#endif
