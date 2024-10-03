#ifndef CIE_FEM_MATHS_CARTESIAN_PRODUCT_HPP
#define CIE_FEM_MATHS_CARTESIAN_PRODUCT_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp" // Ref, Ptr
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <concepts> // std::incrementable


namespace cie::fem::maths {


template <unsigned Dimension>
class OuterProduct
{
public:
    template <concepts::Integer TValue>
    static bool next(unsigned numberOfStates,
                     Ptr<TValue> itStateBegin);

    template <std::incrementable TValue>
    static bool next(Ref<const TValue> firstState,
                     Ref<const TValue> lastState,
                     Ptr<TValue> itStateBegin);

    template <std::incrementable TValue,
              concepts::Iterator<TValue> TItBegin,
              concepts::Iterator<TValue> TItEnd>
    static bool next(TItBegin itOptionsBegin,
                     TItEnd itOptionsEnd,
                     Ptr<TValue> itStateBegin);

    template <std::incrementable TValue>
    static bool next(Ptr<const Ptr<const TValue>> itOptionBeginsArray,
                     Ptr<const Ptr<const TValue>> itOptionEndsArray,
                     Ptr<TValue> itStateBegin);

    constexpr static bool next(Ref<unsigned> rState)
    {return rState++ < (1<<Dimension);}
}; // class OuterProduct


} // namespace cie::fem::maths

#include "packages/maths/impl/OuterProduct_impl.hpp"

#endif
