// --- External Includes ---
#include "benchmark/benchmark.h"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"


template <class TValue, unsigned Dimension>
void linear(benchmark::State& rState)
{
    using Basis = cie::fem::maths::Polynomial<TValue>;
    const cie::fem::maths::AnsatzSpace<Basis,Dimension> ansatzSpace(
        std::vector<Basis> {
            Basis(typename Basis::Coefficients {0.5, -0.5}),
            Basis(typename Basis::Coefficients {0.5,  0.5})
        }
    );

    std::vector<TValue> buffer;
    std::vector<TValue> valueBuffer(ansatzSpace.size());

    using Quadrature = cie::fem::Quadrature<TValue,Dimension>;
    std::vector<typename Quadrature::Point> samples;
    Quadrature(
        cie::fem::GaussLegendreQuadrature<TValue>(/*integrationOrder*/ 5)
    ).getIntegrationPoints(std::back_inserter(samples));

    for ([[maybe_unused]] auto dummy : rState) {
        for (const auto& rSample : samples) {
            ansatzSpace.evaluate(
                rSample,
                valueBuffer);
            benchmark::DoNotOptimize(valueBuffer);
        }
    }
}


template <class TValue, unsigned Dimension>
void quadratic(benchmark::State& rState)
{
    using Basis = cie::fem::maths::Polynomial<TValue>;
    const cie::fem::maths::AnsatzSpace<Basis,Dimension> ansatzSpace(
        std::vector<Basis> {
            Basis(typename Basis::Coefficients {0.5, -0.5}),
            Basis(typename Basis::Coefficients {0.5,  0.5}),
            Basis(typename Basis::Coefficients {1.0,  0.0, -1.0})
        }
    );

    std::vector<TValue> valueBuffer(ansatzSpace.size());

    using Quadrature = cie::fem::Quadrature<TValue,Dimension>;
    std::vector<typename Quadrature::Point> samples;
    Quadrature(
        cie::fem::GaussLegendreQuadrature<TValue>(/*integrationOrder*/ 5)
    ).getIntegrationPoints(std::back_inserter(samples));

    for ([[maybe_unused]] auto dummy : rState) {
        for (const auto& rSample : samples) {
            ansatzSpace.evaluate(rSample.data(),
                                 rSample.data() + rSample.size(),
                                 valueBuffer.data());
            benchmark::DoNotOptimize(valueBuffer);
        }
    }
}


BENCHMARK_TEMPLATE(linear, float, 1)->Name("AnsatzSpace_linear_float_1D");
BENCHMARK_TEMPLATE(linear, float, 2)->Name("AnsatzSpace_linear_float_2D");
BENCHMARK_TEMPLATE(linear, float, 3)->Name("AnsatzSpace_linear_float_3D");

BENCHMARK_TEMPLATE(quadratic, float, 1)->Name("AnsatzSpace_quadratic_float_1D");
BENCHMARK_TEMPLATE(quadratic, float, 2)->Name("AnsatzSpace_quadratic_float_2D");
BENCHMARK_TEMPLATE(quadratic, float, 3)->Name("AnsatzSpace_quadratic_float_3D");

BENCHMARK_TEMPLATE(linear, double, 1)->Name("AnsatzSpace_linear_double_1D");
BENCHMARK_TEMPLATE(linear, double, 2)->Name("AnsatzSpace_linear_double_2D");
BENCHMARK_TEMPLATE(linear, double, 3)->Name("AnsatzSpace_linear_double_3D");

BENCHMARK_TEMPLATE(quadratic, double, 1)->Name("AnsatzSpace_quadratic_double_1D");
BENCHMARK_TEMPLATE(quadratic, double, 2)->Name("AnsatzSpace_quadratic_double_2D");
BENCHMARK_TEMPLATE(quadratic, double, 3)->Name("AnsatzSpace_quadratic_double_3D");
