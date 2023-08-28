#ifndef CIE_UTILS_CHECK_MACROS_HPP
#define CIE_UTILS_CHECK_MACROS_HPP

// --- Internal Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/detail.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <sstream>


/* The preprocessor variables can be defined/undefined through CMake */


#define CIE_CHECK(BOOL_EXPRESSION, MESSAGE) \
    if (!(BOOL_EXPRESSION))  [[unlikely]]   \
    {                                       \
        std::stringstream stream;           \
        stream << MESSAGE;                  \
        CIE_THROW(Exception, stream.str())  \
    }


#define CIE_CHECK_POINTER(pointer)           \
    if (pointer == nullptr) [[unlikely]]     \
    {                                        \
        CIE_THROW(cie::NullPtrException, "") \
    }


#ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
    #define CIE_DETAIL_GET_OUT_OF_RANGE_MACRO(_1, _2, MACRO_NAME, ...) MACRO_NAME

    #define CIE_DETAIL_OUT_OF_RANGE_CHECK_DEFAULT_MESSAGE(BOOL_EXPRESSION) \
        if (!(BOOL_EXPRESSION)) [[unlikely]]                               \
            {CIE_THROW(cie::OutOfRangeException, #BOOL_EXPRESSION)}

    #define CIE_DETAIL_OUT_OF_RANGE_CHECK_WITH_MESSAGE(BOOL_EXPRESSION, MESSAGE) \
        if (!(BOOL_EXPRESSION)) [[unlikely]]                                     \
        {                                                                        \
            std::stringstream stream;                                            \
            stream << MESSAGE;                                                   \
            CIE_THROW(cie::OutOfRangeException, stream.str())                    \
        }

    #define CIE_OUT_OF_RANGE_CHECK(...)                         \
        CIE_EXPAND_ARGUMENTS(                                   \
            CIE_DETAIL_GET_OUT_OF_RANGE_MACRO(                  \
                __VA_ARGS__,                                    \
                CIE_DETAIL_OUT_OF_RANGE_CHECK_WITH_MESSAGE,     \
                CIE_DETAIL_OUT_OF_RANGE_CHECK_DEFAULT_MESSAGE,) \
            (__VA_ARGS__)                                       \
        )

#else
    #define CIE_OUT_OF_RANGE_CHECK(...)
#endif


#ifdef CIE_ENABLE_DIVISION_BY_ZERO_CHECKS
    #define CIE_DIVISION_BY_ZERO_CHECK(BOOL_EXPRESSION)            \
        if (!(BOOL_EXPRESSION)) [[unlikely]]                       \
        {                                                          \
            std::stringstream stream;                              \
            stream << #BOOL_EXPRESSION;                            \
            CIE_THROW(cie::DivisionByZeroException, stream.str()); \
        }
#else
    #define CIE_DIVISION_BY_ZERO_CHECK(BOOL_EXPRESSION) {}
#endif


#ifndef NDEBUG
    #define CIE_DEBUG_CHECK(BOOL_EXPRESSION, MESSAGE)       \
        {                                                   \
            if (!(BOOL_EXPRESSION)) [[unlikely]]            \
            {                                               \
                std::stringstream stream;                   \
                stream << #BOOL_EXPRESSION;                 \
                stream << " (";                             \
                stream << MESSAGE;                          \
                stream << ")";                              \
                CIE_THROW(cie::Exception, stream.str());    \
            }                                               \
        }
#else
    #define CIE_DEBUG_CHECK(BOOL_EXPRESSION, MESSAGE)
#endif


#endif
