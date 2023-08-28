#ifndef CIE_FEM_MATHS_LAGRANGE_POLYNOMIAL_HPP
#define CIE_FEM_MATHS_LAGRANGE_POLYNOMIAL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::fem::maths {


template <class TValue>
class LagrangePolynomial final : public Polynomial<TValue>
{
public:
    LagrangePolynomial() noexcept = default;

    /** @brief Construct a Lagrange polynomial on the provided nodes.
     *  @details The constructed polynomial is of minimal degree that satisfied the following criteria:
     *           - evaluates to 1 at the base node
     *           - vanishes at all other nodes
     *  @param p_nodeBegin: iterator to the first node.
     *  @param p_nodeEnd: iterator past the last node.
     *  @param i_base: index of the base node within the provided node range
     *                 (the polynomial evaluates to 1 at this node and vanishes at the rest).
     */
    LagrangePolynomial(Ptr<const TValue> p_nodeBegin,
                       Ptr<const TValue> p_nodeEnd,
                       Size i_base);
}; // class LagrangePolynomial


} // namespace cie::fem::maths

#endif
