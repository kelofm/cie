// --- FEM Includes ---
#include "packages/maths/inc/LagrangePolynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/NChooseK.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <iterator>


namespace cie::fem::maths {


template <class TValue>
LagrangePolynomial<TValue>::LagrangePolynomial(Ptr<const TValue> p_nodeBegin,
                                               Ptr<const TValue> p_nodeEnd,
                                               Size i_base)
    : Polynomial<TValue>()
{
    CIE_BEGIN_EXCEPTION_TRACING

    const Size polynomialOrder = std::distance(p_nodeBegin, p_nodeEnd);

    if (polynomialOrder) [[likely]] {
        DynamicArray<TValue> localNodes;
        localNodes.reserve(polynomialOrder - 1);
        this->_coefficients.resize(polynomialOrder);

        // Get the base node
        CIE_OUT_OF_RANGE_CHECK(i_base < polynomialOrder)
        const auto it_baseNode = p_nodeBegin + i_base;
        const TValue baseNode = *it_baseNode;

        // Collect and negate the rest.
        // At the same time, compute the denominator
        TValue denominator = static_cast<TValue>(1);
        const auto transform = [&denominator, baseNode](TValue node) {
            denominator *= baseNode - node;
            return -node;
        };
        CIE_DIVISION_BY_ZERO_CHECK(denominator)

        auto it_local = std::back_inserter(localNodes);
        std::transform(p_nodeBegin, it_baseNode, it_local, transform);

        if (i_base < polynomialOrder - 1)
            std::transform(it_baseNode + 1, p_nodeEnd, it_local, transform);

        // Loop backward through the exponents
        const auto it_localNodeBegin = localNodes.begin();
        for (Size numberOfSelectedNodes=0; numberOfSelectedNodes<polynomialOrder; ++numberOfSelectedNodes) {
            const Size exponent = polynomialOrder - 1 - numberOfSelectedNodes;
            utils::NChooseK permutation(polynomialOrder - 1, numberOfSelectedNodes);
            TValue coefficient = 0;
            do {
                TValue term = static_cast<TValue>(1);
                for (auto i_node : permutation)
                    term *= *(it_localNodeBegin + i_node);
                coefficient += term;
            } while (permutation.next());

            this->_coefficients[exponent] = coefficient / denominator;
        } // for numberOfSelectedNodes in range(polynomialOrder)
    } // if polynomialOrder

    CIE_END_EXCEPTION_TRACING
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(LagrangePolynomial);


} // namespace cie::fem::maths
