#ifndef CIE_FEM_MATHS_CARTESIAN_PRODUCT_HPP
#define CIE_FEM_MATHS_CARTESIAN_PRODUCT_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem::maths {


template <unsigned Dimension>
class CartesianProduct
{
public:
    template <concepts::Integer TValue>
    static bool next(unsigned numberOfStates,
                     Ptr<TValue> itStateBegin);

    constexpr static bool next(Ref<unsigned> rState)
    {return rState++ < (1<<Dimension);}
}; // class CartesianProduct


} // namespace cie::fem::maths

#include "packages/maths/impl/CartesianProduct_impl.hpp"

#endif
