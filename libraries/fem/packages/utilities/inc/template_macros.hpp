#ifndef CIE_FEM_TEMPLATE_MACROS_HPP
#define CIE_FEM_TEMPLATE_MACROS_HPP

/**
 * A collection of macros for instantiating templates
 * with specified numeric types and dimensions.
 */

// --- Utility Includes ---
#include "packages/macros/inc/detail.hpp"


#define CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(CLASS_NAME)                    \
    template class CLASS_NAME<float>;                                       \
    template class CLASS_NAME<double>


#define CIE_FEM_INSTANTIATE_TEMPLATE_DIMENSIONS(CLASS_NAME, NUMERIC_TYPE)   \
    template class CLASS_NAME<NUMERIC_TYPE, 1>;                             \
    template class CLASS_NAME<NUMERIC_TYPE, 2>;                             \
    template class CLASS_NAME<NUMERIC_TYPE, 3>;


#define CIE_FEM_INSTANTIATE_TEMPLATE(CLASS_NAME)                            \
    CIE_FEM_INSTANTIATE_TEMPLATE_DIMENSIONS(CLASS_NAME, float)              \
    CIE_FEM_INSTANTIATE_TEMPLATE_DIMENSIONS(CLASS_NAME, double)


#endif