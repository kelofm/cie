// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- FEM Incldues ---
#include "packages/utilities/inc/template_macros.hpp"

// --- Internal Includes ---
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"

// --- STL Includes ---
#include <cmath>
#include <numbers>
#include <functional>
#include <limits>


namespace cie::fem {


template <class T>
struct GaussLegendreInitializer
{
    static std::pair<typename GaussLegendreQuadrature<T>::NodeContainer,
                     typename GaussLegendreQuadrature<T>::WeightContainer>
    getNodesAndWeights(Size integrationOrder,
                       T maxAbsoluteError,
                       Size maxIterations)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        CIE_CHECK(0 < maxAbsoluteError, "the maximum absolute node error must positive")
        CIE_CHECK(0 < maxIterations, "the maximum number of Newton iterations must be at least 1")

        std::pair<typename GaussLegendreQuadrature<T>::NodeContainer, typename GaussLegendreQuadrature<T>::WeightContainer> nodesAndWeights;
        auto& r_nodes   = nodesAndWeights.first;
        auto& r_weights = nodesAndWeights.second;

        r_nodes.resize(integrationOrder);
        r_weights.resize(integrationOrder);

        T legendreValue      = std::numeric_limits<T>::max();
        T legendreDerivative = std::numeric_limits<T>::max();

        // Handle symmetry
        Size i_pivot;
        if (integrationOrder % 2 == 1) {
            i_pivot = integrationOrder / 2 + 1;
            r_nodes[i_pivot] = 0;

            T weight = integrationOrder * std::legendre(integrationOrder-1, T(0));
            weight = 2 / weight / weight;
            r_weights[i_pivot] = weight;

            --i_pivot;
        } else {
            i_pivot = integrationOrder / 2;
        }

        // Compute nodes and weights
        const Size integrationOrderM1 = integrationOrder - 1;

        for (Size index=0; index<=i_pivot; ++index) {
            bool converged = false;
            T node        = approximateLegendreRoot(integrationOrder, index);
            T oneMinusNode2 = 1 - node*node;

            // Newton iteration
            for (Size iteration=0; iteration<maxIterations; ++iteration) {
                // Update helpers
                legendreValue = std::legendre(integrationOrder, node);

                legendreDerivative = integrationOrder * (-node*legendreValue + std::legendre(integrationOrderM1, node));
                legendreDerivative /= oneMinusNode2;

                // Check convergence
                if (std::abs(legendreValue) < maxAbsoluteError) { // TODO: this is incorrect => check node position error instead of value
                    converged = true;
                    break;
                }

                // Update
                node          -= legendreValue / legendreDerivative;
                oneMinusNode2 = 1 - node*node;
            }

            if (!converged) {
                CIE_THROW(Exception, "Computation of Gauss-Legendre nodes failed to converge!")
            }

            T weight = 2 / legendreDerivative / legendreDerivative / oneMinusNode2;

            Size symmetricIndex       = integrationOrder - index - 1;
            r_nodes[index]            = node;
            r_weights[index]          = weight;
            r_nodes[symmetricIndex]   = -node;
            r_weights[symmetricIndex] = weight;
        }

        return nodesAndWeights;

        CIE_END_EXCEPTION_TRACING
    }

    static T approximateLegendreRoot(Size order, Size index)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        index = order - index;
        const T eighth = T(1) / T(8);

        T output = T(1) - (eighth - eighth/order)/(order*order);
        output *= std::cos(std::numbers::pi * T(4*index-1)/T(4*order + 2));

        return output;

        CIE_END_EXCEPTION_TRACING
    }
};


template <concepts::Numeric NT>
GaussLegendreQuadrature<NT>::GaussLegendreQuadrature(Size integrationOrder,
                                                     NT maxAbsoluteNodeError,
                                                     Size maxNewtonIterations)
    : QuadratureBase<NT>(GaussLegendreInitializer<NT>::getNodesAndWeights(integrationOrder,
                                                                          maxAbsoluteNodeError,
                                                                          maxNewtonIterations))
{
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(GaussLegendreQuadrature);


} // namespace cie::fem
