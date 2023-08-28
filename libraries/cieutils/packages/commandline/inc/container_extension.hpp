#ifndef CIE_UTILS_COMMAND_LINE_CONTAINER_EXTENSION_HPP
#define CIE_UTILS_COMMAND_LINE_CONTAINER_EXTENSION_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"


namespace cie::utils {


///@addtogroup cieutils
///@{

// TODO: GCC has an issue with the forward declaration of constrained partial class template specializations
//       => uncomment this and the implementation in the corresponding impl file.
//template <concepts::ReservableContainer TContainer>
//struct ArgParse::Results::ValueConverter<TContainer>;
//
//
//template <concepts::FixedSizeContainer TContainer>
//struct ArgParse::Results::ValueConverter<TContainer>;

template <concepts::Container TContainer>
struct ArgParse::Results::ValueConverter<TContainer>;

///@}


} // namespace cie::utils

#include "packages/commandline/impl/container_extension_impl.hpp"

#endif