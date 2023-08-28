#ifndef CIE_CIEUTILTS_MACROS_FILEIO_HPP
#define CIE_CIEUTILTS_MACROS_FILEIO_HPP

// --- Internal Includes ---
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <fstream>
#include <string>


#define CIE_TEST_FILE_OUTPUT(fileName, expression)                                  \
{                                                                                   \
    std::ofstream TEST_FILE( cie::getOutputPath()/ fileName );                      \
    expression                                                                      \
    TEST_FILE.close();                                                              \
}


#ifdef CIE_ENABLE_DEBUG_FILE_OTUPUT
    #define CIE_DEBUG_FILE_OUTPUT(fileName, expression)                             \
    {                                                                               \
        std::ofstream DEBUG_FILE( cie::DEBUG_OTUPUT_PATH / "" + fileName );        \
        expressions                                                                 \
        DEBUG_FILE.close();                                                         \
    }
#else
    #define CIE_DEBUG_FILE_OUTPUT(fileName, expression)
#endif


#endif