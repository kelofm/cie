#ifndef CIE_UTILS_CONCEPTS_STREAMABLE_HPP
#define CIE_UTILS_CONCEPTS_STREAMABLE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <sstream>
#include <iosfwd>


namespace cie::concepts {


namespace detail {
/** Types that can be directly pushed into string streams.
 *  (excluding function templates that act on streams)
 */
template <class T>
concept StringStreamableType
= AnyOf<typename std::decay<T>::type, short,
                                      unsigned short,
                                      int,
                                      unsigned int,
                                      long,
                                      unsigned long,
                                      long long,
                                      unsigned long long,
                                      float,
                                      double,
                                      long double,
                                      bool,
                                      char*,
                                      const char*,
                                      void*,
                                      volatile void*,
                                      std::nullptr_t,
                                      std::ios_base,
                                      std::streambuf,
                                      std::ostream,
                                      std::string>;

// Function template acting on streams
// (TODO: this just checks for function pointers and std::ostream)
template <class T>
concept StreamInsertOperator
= requires (T t, std::ostream& r_stream)
{
    {t(r_stream)} -> std::same_as<std::ostream&>;
};

//template <class TOperator, class TStream>
//concept StreamInsertOperator
//= requires (TOperator op, TStream& r_stream)
//{
//    {op(r_stream)} -> std::same_as<TStream&>;
//};
} // namespace detail

template <class T>
concept StringStreamable
= detail::StringStreamableType<T> || detail::StreamInsertOperator<T>;


} // namespace cie::concepts


#endif