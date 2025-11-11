#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"


namespace cie::fem::maths {


template <class TValue>
class LegendrePolynomial : public Polynomial<TValue>
{
public:
    LegendrePolynomial() noexcept = default;

    LegendrePolynomial(unsigned index);
}; // class LegendrePolynomial


template <class TValue>
class IntegratedLegendrePolynomial : public Polynomial<TValue>
{
public:
    IntegratedLegendrePolynomial() noexcept = default;

    IntegratedLegendrePolynomial(unsigned index);
}; // class IntegratedLegendrePolynomial


} // namespace cie::fem::maths
