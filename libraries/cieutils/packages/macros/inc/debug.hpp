#ifndef NDEBUG


// --- Utility Includes ---
#include "packages/macros/inc/detail.hpp"

// --- STL Includes ---
#include <iostream>


#define CIE_DEBUG_LINE(...) \
    std::cout << __VA_ARGS__ << CIE_CODE_LOCATION << std::endl;

#else


#define CIE_DEBUG_LINE(...)


#endif
