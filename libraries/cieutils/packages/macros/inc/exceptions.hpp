#ifndef CIE_UTILS_EXCEPTION_MACROS_HPP
#define CIE_UTILS_EXCEPTION_MACROS_HPP

// --- Internal Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/detail.hpp"

// --- STL Includes ---
#include <sstream>
#include <concepts>


/* Define exception-related macros */

#define CIE_THROW(TEXCEPTION, MESSAGE)                             \
    {                                                              \
        std::stringstream STREAM;                                  \
        STREAM << MESSAGE;                                         \
        throw cie::exceptionFactory<TEXCEPTION>(CIE_CODE_LOCATION, \
                                                STREAM.str(),      \
                                                #TEXCEPTION);      \
    }


#define CIE_RETHROW(EXCEPTION, ADDITIONAL_MESSAGE)                                   \
    {                                                                                \
        std::stringstream EXCEPTION_STREAM;                                          \
        EXCEPTION_STREAM << ADDITIONAL_MESSAGE;                                      \
        throw cie::exceptionFactory<typename std::decay<decltype(EXCEPTION)>::type>( \
            EXCEPTION, CIE_CODE_LOCATION, EXCEPTION_STREAM.str()                     \
        );                                                                           \
    }


#ifdef CIE_ENABLE_EXCEPTION_TRACING

    #define CIE_BEGIN_EXCEPTION_TRACING try {

    #define CIE_END_EXCEPTION_TRACING                      \
        }                                                  \
        catch (const cie::Exception& EXCEPTION)            \
        {                                                  \
            throw cie::Exception(                          \
                CIE_CODE_LOCATION,                         \
                EXCEPTION.what(),                          \
                EXCEPTION.stackLevel() + 1                 \
            );                                             \
        }                                                  \
        catch (const std::exception& EXCEPTION)            \
        {                                                  \
            throw cie::exceptionFactory<cie::Exception>(   \
                CIE_CODE_LOCATION,                         \
                EXCEPTION.what(),                          \
                "std::exception"                           \
            );                                             \
        }                                                  \
        catch (...)                                        \
        {                                                  \
            CIE_THROW(cie::Exception, "Unknown exception") \
        }

#else

    #define CIE_BEGIN_EXCEPTION_TRACING

    #define CIE_END_EXCEPTION_TRACING

#endif // CIE_ENABLE_EXCEPTION_TRACING


#endif