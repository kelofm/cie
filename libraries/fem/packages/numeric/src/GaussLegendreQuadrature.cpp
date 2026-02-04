#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

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

#if __STDCPP_MATH_SPEC_FUNCS__ < 201003L // I hate clang

#ifndef __clang__
static_assert(false, "terribleLegendre is supposed to replace std::legendre exclusively on clang");
#endif

#include "tsl/robin_map.h"

/** @brief Naive recursive evaluation of a Legendre polynomial.
 *  @details This function directly uses Legendre polynomials'
 *           recurrence relation to compute their values at the
 *           given argument in [-1,1]. The recurrence relation
 *           reads:
 *           @f[
 *              p_0(x) &= 1 \
 *              p_1(x) &= x \
 *              p_k(x) &= \frac{2(k-1)+1}{k} x p_{k-1}(x)
 *                        -  \frac{k-1}{k}p_{k-2}(x)
 *           @f]
 *  @note Both performance and accuracy are terrible. Ideally,
 *        one would use the standard implementation instead, but
 *        clang still hasn't implemented mathematical special
 *        functions that have been part of the standard since C++17.
 *  @todo Replace @a terribleLegendre with @ref std::legendre once
 *        it gets implemented in clang's standard library, or use
 *        the boost implementation instead.
 */
long double terribleLegendre(std::size_t order,
                             long double x,
                             tsl::robin_map<std::size_t,long double>& rCache)
{
    if (order == 0) {
        return 1.0;
    } else if (order == 1) {
        return x;
    } else {
        const auto it = rCache.find(order);
        if (it == rCache.end()) {
            const long double prev = terribleLegendre(order - 1, x, rCache);
            const long double prevPrev = terribleLegendre(order - 2, x, rCache);
            const long double current = (2 * (order - 1) + 1) / (long double)order * x * prev
                                      - (order - 1) / (long double)order * prevPrev;
            rCache.emplace(order, current);
            return current;
        } else {
            return it->second;
        }
    }
}

namespace std {
double legendre(size_t order, double x)
{
    tsl::robin_map<std::size_t,long double> cache;
    return terribleLegendre(order, x, cache);
}
} // namespace std

#endif


namespace cie::fem {


template <class T>
struct GaussLegendreInitializer {
    static std::pair<typename GaussLegendreQuadrature<T>::NodeContainer,
                     typename GaussLegendreQuadrature<T>::WeightContainer>
    getNodesAndWeights(Size integrationOrder,
                       utils::Comparison<T> comparison,
                       Size maxIterations)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        CIE_CHECK(0 < integrationOrder, "the integration order must be positive")
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
        if (integrationOrder % 2 == 1 && 1 < integrationOrder) {
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
                if (comparison.equal(std::abs(legendreValue), T(0))) { /// @todo this is incorrect => check node position error instead of value
                    converged = true;
                    break;
                }

                // Update
                node          -= legendreValue / legendreDerivative;
                oneMinusNode2 = 1 - node*node;
            }

            CIE_CHECK(converged,
                      "Calculation of Gauss-Legendre nodes failed to converge within "
                      << maxIterations << " iterations for node " << index
                      << "/" << integrationOrder
                      << " (" << legendreValue << " residual).\n")

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

    static T approximateLegendreRoot(Size order, Size index) {
        index = order - index;
        const T eighth = T(1) / T(8);

        T output = T(1) - (eighth - eighth/order)/(order*order);
        output *= std::cos(std::numbers::pi * T(4*index-1)/T(4*order + 2));

        return output;
    }
};


template <concepts::Numeric NT>
GaussLegendreQuadrature<NT>::GaussLegendreQuadrature(Size integrationOrder,
                                                     utils::Comparison<NT> comparison,
                                                     Size maxNewtonIterations)
    : QuadratureBase<NT>(GaussLegendreInitializer<NT>::getNodesAndWeights(
        integrationOrder,
        comparison,
        maxNewtonIterations))
{
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(GaussLegendreQuadrature);


} // namespace cie::fem
