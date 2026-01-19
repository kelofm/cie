#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
PolynomialView<TValue>::PolynomialView(ConstSpan coefficients) noexcept
requires (!hasStaticBasis)
    : _coefficients(coefficients)
{}


template <concepts::Numeric TValue>
constexpr PolynomialView<TValue>::PolynomialView(std::span<const TValue,coefficientCount> coefficients) noexcept
requires hasStaticBasis
    : _coefficients(coefficients)
{}


template <concepts::Numeric TValue>
typename PolynomialView<TValue>::Derivative
PolynomialView<TValue>::makeDerivative(Span buffer) const {
    if (!_coefficients.empty() && buffer.size() != (_coefficients.empty() ? 0ul : _coefficients.size() - 1)) {
        CIE_THROW(OutOfRangeException,
                  "required buffer size is " << _coefficients.size() - 1 << " "
                    << "but got " << buffer.size());
    }

    const auto polynomialOrder = _coefficients.size();

    if (1 < polynomialOrder) [[likely]] {
        // Push first coefficient (no multiplication required)
        buffer.front() = _coefficients[1];
        if (2 < polynomialOrder) {
            const auto itCoefficientEnd = _coefficients.end();
            TValue power = static_cast<TValue>(2);
            auto itBuffer = buffer.begin() + 1;
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power, ++itBuffer)
                *itBuffer = power * (*itCoefficient);
        } // if 2 < polynoialOrder
    } // if 1 < polynomialOrder

    return PolynomialView(buffer);
}


template <concepts::Numeric TValue>
void PolynomialView<TValue>::evaluate(ConstSpan in, Span out) const {
    out.front() = utils::evaluatePolynomialHorner(
        in.front(),
        _coefficients.begin(),
        _coefficients.end());
}


template <concepts::Numeric TValue>
unsigned PolynomialView<TValue>::size() const noexcept
{
    return 1u;
}


template <concepts::Numeric TValue>
typename PolynomialView<TValue>::ConstSpan
PolynomialView<TValue>::coefficients() const noexcept {
    return _coefficients;
}


template <class TValue>
Polynomial<TValue>::Polynomial(const Polynomial& rRight)
    : _coefficients(rRight._coefficients),
      _wrapped()
{
    _wrapped = PolynomialView<TValue>(_coefficients);
}


template <class TValue>
Polynomial<TValue>::Polynomial(RightRef<Coefficients> rCoefficients) noexcept
    : _coefficients(std::move(rCoefficients))
{
    _wrapped = PolynomialView<TValue>(_coefficients);
}


template <class TValue>
Polynomial<TValue>::Polynomial(ConstSpan coefficients)
    : _coefficients(coefficients.begin(), coefficients.end())
{
    _wrapped = PolynomialView<TValue>(_coefficients);
}


template <class TValue>
Polynomial<TValue>& Polynomial<TValue>::operator=(const Polynomial& rRight) {
    _coefficients = rRight._coefficients;
    _wrapped = PolynomialView<TValue>(_coefficients);
    return *this;
}


template <class TValue>
void Polynomial<TValue>::evaluate(ConstSpan in, Span out) const
{
    _wrapped.evaluate(in, out);
}


template <class TValue>
unsigned Polynomial<TValue>::size() const noexcept
{
    return _wrapped.size();
}


template <class TValue>
Polynomial<TValue> Polynomial<TValue>::makeDerivative() const {
    Polynomial derivative;
    derivative._coefficients.resize(_coefficients.empty() ? 0 : _coefficients.size() - 1);
    derivative._wrapped = _wrapped.makeDerivative(derivative._coefficients);
    return derivative;
}


template <class TValue>
std::span<const TValue> Polynomial<TValue>::coefficients() const noexcept {
    return {_coefficients.data(), _coefficients.size()};
}


} // namespace cie::fem::maths
