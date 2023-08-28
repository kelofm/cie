#ifndef CIE_UTILS_MACROS_WARNING_HPP
#define CIE_UTILS_MACROS_WARNING_HPP

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"


#define CIE_WARNING( message )                      \
    cie::utils::LoggerSingleton::get().warn( message );


#endif