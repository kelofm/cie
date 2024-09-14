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
     *  @details The constructed polynomial is of minimal degree that satisfies the following criteria:
     *           - evaluates to 1 at the base node
     *           - vanishes at all other nodes
     *  @param pNodeBegin: iterator to the first node.
     *  @param pNodeEnd: iterator past the last node.
     *  @param iBase: index of the base node within the provided node range
     *                 (the polynomial evaluates to 1 at this node and vanishes at the rest).
     */
    LagrangePolynomial(Ptr<const TValue> pNodeBegin,
                       Ptr<const TValue> pNodeEnd,
                       Size iBase);
}; // class LagrangePolynomial


} // namespace cie::fem::maths

#endif
