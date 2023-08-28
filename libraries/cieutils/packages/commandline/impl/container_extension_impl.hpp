#ifndef CIE_UTILS_COMMAND_LINE_CONTAINER_EXTENSION_IMPL_HPP
#define CIE_UTILS_COMMAND_LINE_CONTAINER_EXTENSION_IMPL_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::utils {


//template <concepts::ReservableContainer TContainer>
//struct ArgParse::Results::ValueConverter<TContainer>
//{
//    static TContainer convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
//    {
//        CIE_BEGIN_EXCEPTION_TRACING
//
//        using TValue = typename std::decay<typename TContainer::value_type>::type;
//
//        TContainer output;
//        utils::reserve(output, std::distance(begin, end));
//
//        for (; begin!=end; ++begin)
//        {
//            auto localEnd = begin;
//            ++localEnd;
//            output.emplace_back(ArgParse::Results::ValueConverter<TValue>(begin, localEnd));
//        }
//
//        return output;
//
//        CIE_END_EXCEPTION_TRACING
//    }
//}; // struct ArgParse::Results::ValueConverter<ReservableContainer>
//
//
//template <concepts::FixedSizeContainer TContainer>
//struct ArgParse::Results::ValueConverter<TContainer>
//{
//    static TContainer convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
//    {
//        CIE_BEGIN_EXCEPTION_TRACING
//
//        using TValue = typename std::decay<typename TContainer::value_type>::type;
//
//        TContainer output;
//        CIE_OUT_OF_RANGE_CHECK(output.size() == std::distance(begin, end));
//
//        auto it_output = output.begin();
//        for (; begin!=end; ++begin, ++it_output)
//        {
//            auto localEnd = begin;
//            ++localEnd;
//            *it_output = ArgParse::Results::ValueConverter<TValue>(begin, localEnd);
//        }
//
//        return output;
//
//        CIE_END_EXCEPTION_TRACING
//    }
//}; // struct ArgParse::Results::ValueConverter<ReservableContainer>


namespace detail {

template <concepts::ReservableContainer TContainer>
TContainer convertArgParseResults(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
{
    CIE_BEGIN_EXCEPTION_TRACING

        using TValue = typename std::decay<typename TContainer::value_type>::type;

        TContainer output;
        utils::reserve(output, std::distance(begin, end));

        for (; begin!=end; ++begin)
        {
            auto localEnd = begin;
            ++localEnd;
            output.emplace_back(ArgParse::Results::ValueConverter<TValue>::convert(begin, localEnd));
        }

        return output;

        CIE_END_EXCEPTION_TRACING
};


template <concepts::FixedSizeContainer TContainer>
TContainer convertArgParseResults(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
{
    CIE_BEGIN_EXCEPTION_TRACING

        using TValue = typename std::decay<typename TContainer::value_type>::type;

        TContainer output;
        CIE_OUT_OF_RANGE_CHECK(output.size() == std::distance(begin, end));

        for (; begin!=end; ++begin)
        {
            auto localEnd = begin;
            ++localEnd;
            output.emplace_back(ArgParse::Results::ValueConverter<TValue>(begin, localEnd));
        }

        return output;

        CIE_END_EXCEPTION_TRACING
};

} // namespace detail


template <concepts::Container TContainer>
struct ArgParse::Results::ValueConverter<TContainer>
{
    static TContainer convert(ArgParse::Results::ValueIterator begin, ArgParse::Results::ValueIterator end)
    {
        return detail::convertArgParseResults<TContainer>(begin, end);
    }
};


} // namespace cie::utils


#endif