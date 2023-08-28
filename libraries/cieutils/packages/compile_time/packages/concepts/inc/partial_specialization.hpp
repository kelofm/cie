#ifndef CIE_UTILS_CONCEPTS_PARTIAL_PSECIALIZATION_HPP
#define CIE_UTILS_CONCEPTS_PARTIAL_PSECIALIZATION_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::concepts::detail {


template <class Tag, class ...Args>
requires Tag::value
struct TrueTemplateBase {};


} // namespace cie::concepts::detail


/// Define an invalid class template to be (partially) specialized later 
#define CIE_DEFINE_INVALID_CLASS_TEMPLATE_TO_SPECIALIZE( className )           \
    template <class ...Args>                                                        \
    class className :                                                               \
        public cie::concepts::detail::TrueTemplateBase<std::false_type,Args...>     \
    {};


#endif