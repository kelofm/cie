#ifndef CIE_UTILS_MACROS_TESTING_HPP
#define CIE_UTILS_MACROS_TESTING_HPP

// --- Internal Includes ---
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <string>


#ifndef CIE_TESTRUNNER_NAME
#define CIE_TESTRUNNER_NAME unnamed_testrunner
#endif


#define CIE_DETAIL_STRINGIFY(argument) #argument

#define CIE_STRINGIFY(argument) CIE_DETAIL_STRINGIFY(argument)


#define CIE_TEST_CASE_INIT(testCaseName)    \
    cie::utils::LogBlock LOG_BLOCK(         \
        testCaseName,                       \
        cie::utils::LoggerSingleton::get()  \
    );


#endif