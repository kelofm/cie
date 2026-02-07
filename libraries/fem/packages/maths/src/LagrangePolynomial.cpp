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
LagrangePolynomial<TValue>::LagrangePolynomial(std::span<const TValue> nodes,
                                               Size iBase)
    : Polynomial<TValue>()
{
    CIE_BEGIN_EXCEPTION_TRACING

    const Size polynomialOrder = nodes.size();
    typename LagrangePolynomial::Coefficients coefficients(polynomialOrder);

    if (polynomialOrder) [[likely]] {
        DynamicArray<TValue> localNodes;
        localNodes.reserve(polynomialOrder - 1);

        // Get the base node
        CIE_OUT_OF_RANGE_CHECK(iBase < polynomialOrder)
        const auto itBaseNode = nodes.begin() + iBase;
        const TValue baseNode = *itBaseNode;

        // Collect and negate the rest.
        // At the same time, compute the denominator
        TValue denominator = static_cast<TValue>(1);
        const auto transform = [&denominator, baseNode](TValue node) {
            denominator *= baseNode - node;
            return -node;
        };
        CIE_DIVISION_BY_ZERO_CHECK(denominator)

        auto itLocal = std::back_inserter(localNodes);
        std::transform(nodes.begin(), itBaseNode, itLocal, transform);

        if (iBase < polynomialOrder - 1)
            std::transform(itBaseNode + 1, nodes.end(), itLocal, transform);

        // Loop backward through the exponents
        const auto itLocalNodeBegin = localNodes.begin();
        for (Size numberOfSelectedNodes=0; numberOfSelectedNodes<polynomialOrder; ++numberOfSelectedNodes) {
            const Size exponent = polynomialOrder - 1 - numberOfSelectedNodes;
            utils::NChooseK permutation(polynomialOrder - 1, numberOfSelectedNodes);
            TValue coefficient = 0;
            do {
                TValue term = static_cast<TValue>(1);
                for (auto iNode : permutation)
                    term *= *(itLocalNodeBegin + iNode);
                coefficient += term;
            } while (permutation.next());

            coefficients[exponent] = coefficient / denominator;
        } // for numberOfSelectedNodes in range(polynomialOrder)
    } // if polynomialOrder

    static_cast<Polynomial<TValue>&>(*this) = Polynomial<TValue>(std::move(coefficients));

    CIE_END_EXCEPTION_TRACING
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(LagrangePolynomial);


} // namespace cie::fem::maths
