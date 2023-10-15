#ifndef CIE_UTILS_TYPE_SELECT_INTEGER_HPP
#define CIE_UTILS_TYPE_SELECT_INTEGER_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Concatenate.hpp"

// --- STL Includes ---
#include <limits>
#include <tuple>
#include <cstdint>


namespace cie::ct {


namespace detail {
template <class TIndex>
class IntegerSelector
{
public:
    /// @brief Select the shortest integer type from the provided tuple that can represent the given range.
    ///
    /// @tparam Min: range begin.
    /// @tparam Max: range end.
    /// @tparam TTypeTuple: tuple with integer types sorted in ascending order by size.
    /// @note @a Type is @a void if none of the integer types in @p TTypeTuple can represent the specified range.
    template <TIndex Min, TIndex Max, class TTypeTuple>
    class ShortestWithRange
    {
    private:
        using ExtendedTuple = typename ct::ConcatenateTuple<TTypeTuple>::template Type<std::tuple<void>>;

        template <Size I>
        requires (I < std::tuple_size_v<ExtendedTuple> - 1)
        static constexpr Size getIndex()
        {
            using Current = std::tuple_element_t<I,ExtendedTuple>;
            constexpr const Current min = std::numeric_limits<Current>::min();
            constexpr const Current max = std::numeric_limits<Current>::max();
            if constexpr (min <= Min && Max <= max)
                return I;
            else
                return getIndex<I+1>();
        }

        template <Size I>
        requires (I == std::tuple_size_v<ExtendedTuple> - 1)
        static constexpr Size getIndex()
        {
            return std::tuple_size_v<ExtendedTuple> - 1;
        }

    public:
        using Type = std::tuple_element_t<getIndex<0>(),ExtendedTuple>;
    }; // class ShortestWithRange
}; // class IntegerSelector

using SignedIntegers = std::tuple<int8_t,
                                  int16_t,
                                  int32_t,
                                  int64_t,
                                  int_fast8_t,
                                  int_fast16_t,
                                  int_fast32_t,
                                  int_fast64_t,
                                  int_least8_t,
                                  int_least16_t,
                                  int_least32_t,
                                  int_least64_t>;

using UnsignedIntegers = std::tuple<uint_least8_t,
                                    uint8_t,
                                    uint_fast8_t,
                                    uint_least16_t,
                                    uint16_t,
                                    uint_fast16_t,
                                    uint_least32_t,
                                    uint32_t,
                                    uint_fast32_t,
                                    uint_least64_t,
                                    uint64_t,
                                    uint_fast64_t>;
} // namespace detail


/// @addtogroup cieutils
struct Signed
{
    using Longest = intmax_t;

    template <Longest Min, Longest Max>
    using ShortestWithRange = typename detail::IntegerSelector<Longest>::ShortestWithRange<Min,Max,detail::SignedIntegers>::Type;


}; // struct SelectSigned


/// @addtogroup cieutils
struct Unsigned
{
    using Longest = uintmax_t;

    template <Longest Min, Longest Max>
    using ShortestWithRange = typename detail::IntegerSelector<Longest>::ShortestWithRange<Min,Max,detail::UnsignedIntegers>::Type;
}; // struct SelectUnsigned


} // namespace cie::ct


#endif
