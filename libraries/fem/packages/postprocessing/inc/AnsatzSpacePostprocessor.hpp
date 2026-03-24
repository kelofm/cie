#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"

// --- STL Includes ---
#include <span>


namespace cie::fem {


template <class TValue, unsigned Dimension, std::size_t SetSize = 0ul>
class AnsatzSpacePostprocessor {
public:
    using Ansatz = maths::AnsatzSpace<
        maths::Polynomial<TValue>,
        Dimension>;

    //static void interpolate(
    //    Ref<const >)
}; // class AnsatzSpacePostprocessor


} // namespace cie::fem
