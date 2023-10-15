#ifndef CIE_UTILS_TYPES_MODIFIERS_HPP
#define CIE_UTILS_TYPES_MODIFIERS_HPP

// --- Utility Inlcludes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <type_traits>


namespace cie {


/// @addtogroup cieutils
template <class Source, class Target>
struct CopyConstQualifier
{using Type = Target;};


/// @addtogroup cieutils
template <concepts::Const Source, class Target>
struct CopyConstQualifier<Source, Target>
{using Type = const Target;};


} // namespace cie


#endif