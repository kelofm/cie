// --- External Includes ---
#include "benchmark/benchmark.h"

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"


template <class TValue, unsigned Dimension>
void linear(benchmark::State& r_state)
{
    using Basis = cie::fem::maths::Polynomial<TValue>;
    const cie::fem::maths::AnsatzSpace<Basis,Dimension> ansatzSpace(
        cie::DynamicArray<Basis> {
            Basis(typename Basis::Coefficients {0.5, -0.5}),
            Basis(typename Basis::Coefficients {0.5,  0.5})
        }
    );

    cie::DynamicArray<TValue> valueBuffer(ansatzSpace.size());

    using Quadrature = cie::fem::Quadrature<TValue,Dimension>;
    cie::DynamicArray<typename Quadrature::Point> samples;
    Quadrature(
        cie::fem::GaussLegendreQuadrature<TValue>(/*integrationOrder*/ 5)
    ).getIntegrationPoints(std::back_inserter(samples));

    for ([[maybe_unused]] auto dummy : r_state) {
        for (const auto& r_sample : samples) {
            ansatzSpace.evaluate(r_sample.data(),
                                 r_sample.data() + r_sample.size(),
                                 valueBuffer.data());
            benchmark::DoNotOptimize(valueBuffer);
        }
    }
}


template <class TValue, unsigned Dimension>
void quadratic(benchmark::State& r_state)
{
    using Basis = cie::fem::maths::Polynomial<TValue>;
    const cie::fem::maths::AnsatzSpace<Basis,Dimension> ansatzSpace(
        cie::DynamicArray<Basis> {
            Basis(typename Basis::Coefficients {0.5, -0.5}),
            Basis(typename Basis::Coefficients {0.5,  0.5}),
            Basis(typename Basis::Coefficients {1.0,  0.0, -1.0})
        }
    );

    cie::DynamicArray<TValue> valueBuffer(ansatzSpace.size());

    using Quadrature = cie::fem::Quadrature<TValue,Dimension>;
    cie::DynamicArray<typename Quadrature::Point> samples;
    Quadrature(
        cie::fem::GaussLegendreQuadrature<TValue>(/*integrationOrder*/ 5)
    ).getIntegrationPoints(std::back_inserter(samples));

    for ([[maybe_unused]] auto dummy : r_state) {
        for (const auto& r_sample : samples) {
            ansatzSpace.evaluate(r_sample.data(),
                                 r_sample.data() + r_sample.size(),
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
