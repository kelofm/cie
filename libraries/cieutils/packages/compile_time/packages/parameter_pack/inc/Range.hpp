#ifndef CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_RANGE_HPP
#define CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_RANGE_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"

// --- STL Includes ---
#include <type_traits>
#include <tuple>
#include <utility>


namespace cie::ct {


///@addtogroup cieutils
///@{

template <class ...Ts>
class Range
{
private:
    using Tuple = std::tuple<Ts...>;

public:
    template <Size Begin, Size End>
    class Subrange {};

    template <Size Begin, Size End>
    requires (Begin <= End && Begin < sizeof...(Ts) && End <= sizeof...(Ts))
    class Subrange<Begin,End>
    {
    private:
        using IndexRange = std::make_integer_sequence<Size,End-Begin>;

        template <class T>
        struct TupleFactory {};

        template <Size ...Is>
        struct TupleFactory<std::integer_sequence<Size,Is...>>
        {using Type = std::tuple<typename std::tuple_element<Is+Begin,Range::Tuple>::type...>;};

        template <template <class ...> class TT, class ...TTs>
        struct ApplyFactory {};

        template <template <class ...> class TT, class ...TTs>
        struct ApplyFactory<TT, std::tuple<TTs...>>
        {using Type = TT<TTs...>;};

    public:
        /// Tuple containing the specified subrange of types
        using Tuple = typename TupleFactory<IndexRange>::Type;

        /// Insert the subrange of types into a nested template
        template <template <class ...> class TT, class ...TArgs>
        using Apply = typename ApplyFactory<TT,Tuple>::Type;
    }; // class Subrange

    template <Size Begin, Size End>
    requires (Begin == End && End == PackSize<Ts...>)
    class Subrange<Begin,End>
    {
    public:
        using Tuple = typename std::tuple<>;

        template <template <class ...> class TT, class ...TArgs>
        using ApplyFactory = TT<>;
    }; // class Subrange<Begin==End>

    /// Insert the range of types into a nested template
    template <template <class ...> class TT, class ...TArgs>
    using Apply = typename Subrange<0,PackSize<Ts...>-1>::template Apply<TT,TArgs...>;
}; // class Range

///@}

} // namespace cie::ct


#endif