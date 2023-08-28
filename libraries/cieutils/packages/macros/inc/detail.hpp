#ifndef CIE_CIEUTILS_MACRO_DETAILS_HPP
#define CIE_CIEUTILS_MACRO_DETAILS_HPP


// File path and line number
#define CIE_MACRO_TO_STRING_HELPER(num) #num
#define CIE_MACRO_TO_STRING(num) CIE_MACRO_TO_STRING_HELPER(num)
#define CIE_CODE_LOCATION __FILE__ ":" CIE_MACRO_TO_STRING(__LINE__)

// Argument unpacking fix for the buggy MSVC preprocessor
#define CIE_EXPAND_ARGUMENTS(arguments) arguments


#endif
