#ifndef CIE_UTILS_COMPILE_TIME_SEQUENCE_HPP
#define CIE_UTILS_COMPILE_TIME_SEQUENCE_HPP

// --- STL Includes ---
#include <utility>


namespace cie::ct {


namespace impl {
template <std::size_t, class>
struct ReverseIndexSequence;

template <std::size_t N, std::size_t ...Is>
struct ReverseIndexSequence<N, std::index_sequence<Is...>>
    : public std::index_sequence<(N - Is)...>
{
};
} // namespace impl


/// @ingroup cieutils
template <std::size_t N>
using MakeIndexSequence = std::make_index_sequence<N>;


/// @ingroup cieutils
template <std::size_t ...Is>
using IndexSequence = std::index_sequence<Is...>;


/// @ingroup cieutils
template <std::size_t N>
requires (0 < N)
using MakeReverseIndexSequence = impl::ReverseIndexSequence<N-1, std::make_index_sequence<N>>;


/// @ingroup cieutils
template <std::size_t ...Is>
using ReverseIndexSequence = IndexSequence<Is...>;


} // namespace cie::ct


#endif
