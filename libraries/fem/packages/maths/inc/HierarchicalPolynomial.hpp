#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::fem::maths {


template <class TValue>
class HierarchicalPolynomial : public Polynomial<TValue>
{
public:
    HierarchicalPolynomial() noexcept = default;

    HierarchicalPolynomial(unsigned index);
}; // class HierarchicalPolynomial


} // namespace cie::fem::maths
