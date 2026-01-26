#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr PolynomialView<TValue,PolynomialOrder>::PolynomialView() noexcept
    : _coefficients(static_cast<const TValue*>(nullptr), coefficientCount)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
PolynomialView<TValue,PolynomialOrder>::PolynomialView(ConstSpan coefficients) noexcept
requires (!hasStaticCoefficients)
    : _coefficients(coefficients)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr PolynomialView<TValue,PolynomialOrder>::PolynomialView(std::span<const TValue,coefficientCount> coefficients) noexcept
requires (hasStaticCoefficients)
    : _coefficients(coefficients)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
typename PolynomialView<TValue,PolynomialOrder>::Derivative
PolynomialView<TValue,PolynomialOrder>::makeDerivative(Span buffer) const
requires (!hasStaticCoefficients) {
    if (!_coefficients.empty() && buffer.size() != (_coefficients.empty() ? 0ul : _coefficients.size() - 1)) {
        CIE_THROW(OutOfRangeException,
                  "required buffer size is " << _coefficients.size() - 1 << " "
                    << "but got " << buffer.size());
    }

    const auto coefficientCount = _coefficients.size();

    if (1 < coefficientCount) [[likely]] {
        // Push first coefficient (no multiplication required)
        buffer.front() = _coefficients[1];
        if (2 < coefficientCount) {
            const auto itCoefficientEnd = _coefficients.end();
            TValue power = static_cast<TValue>(2);
            auto itBuffer = buffer.begin() + 1;
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power, ++itBuffer)
                *itBuffer = power * (*itCoefficient);
        } // if 2 < polynoialOrder
    } // if 1 < coefficientCount

    return PolynomialView(buffer);
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr typename PolynomialView<TValue,PolynomialOrder>::Derivative
PolynomialView<TValue,PolynomialOrder>::makeDerivative(std::span<TValue,Derivative::coefficientCount> buffer) const noexcept
requires (hasStaticCoefficients) {
    if constexpr (1 < coefficientCount) {
        // Push first coefficient (no multiplication required)
        buffer.front() = _coefficients[1];
        if constexpr (2 < coefficientCount) {
            const auto itCoefficientEnd = _coefficients.end();
            TValue power = static_cast<TValue>(2);
            auto itBuffer = buffer.begin() + 1;
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power, ++itBuffer)
                *itBuffer = power * (*itCoefficient);
        } // if 2 < polynoialOrder
    } /* if 1 < PolynomialOrder */

    return Derivative(buffer);
}


template <concepts::Numeric TValue, int PolynomialOrder>
void PolynomialView<TValue,PolynomialOrder>::evaluate(ConstSpan in, Span out) const
requires (!hasStaticCoefficients) {
    out.front() = utils::evaluatePolynomial<utils::PolynomialEvaluation::Horner>(
        in.front(),
        _coefficients);
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr void PolynomialView<TValue,PolynomialOrder>::evaluate(ConstSpan in, Span out) const
requires (hasStaticCoefficients) {
    out.front() = utils::evaluatePolynomial<utils::PolynomialEvaluation::Horner>(
        in.front(),
        _coefficients);
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr unsigned PolynomialView<TValue,PolynomialOrder>::size() noexcept {
    return 1u;
}


template <concepts::Numeric TValue, int PolynomialOrder>
typename PolynomialView<TValue,PolynomialOrder>::ConstSpan
PolynomialView<TValue,PolynomialOrder>::coefficients() const noexcept
requires (!hasStaticCoefficients) {
    return _coefficients;
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr std::span<const TValue,PolynomialView<TValue,PolynomialOrder>::coefficientCount>
PolynomialView<TValue,PolynomialOrder>::coefficients() const noexcept
requires (hasStaticCoefficients) {
    return _coefficients;
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>::Polynomial(RightRef<Coefficients> rCoefficients) noexcept
    : _coefficients(std::move(rCoefficients))
{}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>::Polynomial(Polynomial&& rRight) noexcept
    : _coefficients(std::move(rRight._coefficients))
{}


template <concepts::Numeric TValue, int PolynomialOrder>
Polynomial<TValue,PolynomialOrder>::Polynomial(Polynomial&& rRight)
requires (!hasStaticCoefficients)
    : _coefficients(std::move(rRight._coefficients))
{}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>::Polynomial(Polynomial&& rRight)
requires (hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
Polynomial<TValue,PolynomialOrder>::Polynomial(const Polynomial& rRight)
requires (!hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>::Polynomial(const Polynomial& rRight)
requires (hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric TValue, int PolynomialOrder>
Polynomial<TValue,PolynomialOrder>::Polynomial(ConstSpan coefficients)
requires (!hasStaticCoefficients)
    : _coefficients(coefficients.begin(), coefficients.end())
{}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>::Polynomial(std::span<const TValue,coefficientCount> coefficients)
requires (hasStaticCoefficients)
    : _coefficients()
{
    std::copy_n(
        coefficients.data(),
        coefficientCount,
        _coefficients.data());
}


template <concepts::Numeric TValue, int PolynomialOrder>
Polynomial<TValue,PolynomialOrder>&
Polynomial<TValue,PolynomialOrder>::operator=(Polynomial&& rRight) noexcept
requires (!hasStaticCoefficients) {
    _coefficients = std::move(rRight._coefficients);
    return *this;
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>&
Polynomial<TValue,PolynomialOrder>::operator=(Polynomial&& rRight) noexcept
requires (hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric TValue, int PolynomialOrder>
Polynomial<TValue,PolynomialOrder>& Polynomial<TValue,PolynomialOrder>::operator=(const Polynomial& rRight)
requires (!hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr Polynomial<TValue,PolynomialOrder>& Polynomial<TValue,PolynomialOrder>::operator=(const Polynomial& rRight) noexcept
requires (hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric TValue, int PolynomialOrder>
void Polynomial<TValue,PolynomialOrder>::evaluate(ConstSpan in, Span out) const
requires (!hasStaticCoefficients) {
    this->makeView().evaluate(in, out);
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr void Polynomial<TValue,PolynomialOrder>::evaluate(ConstSpan in, Span out) const noexcept
requires (hasStaticCoefficients) {
    this->makeView().evaluate(in, out);
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr unsigned Polynomial<TValue,PolynomialOrder>::size() noexcept {
    return View::size();
}


template <concepts::Numeric TValue, int PolynomialOrder>
typename Polynomial<TValue,PolynomialOrder>::Derivative
Polynomial<TValue,PolynomialOrder>::makeDerivative() const
requires (!hasStaticCoefficients) {
    Polynomial derivative;
    derivative._coefficients.resize(_coefficients.empty() ? 0 : _coefficients.size() - 1);
    this->makeView().makeDerivative(derivative._coefficients);
    return derivative;
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr typename Polynomial<TValue,PolynomialOrder>::Derivative
Polynomial<TValue,PolynomialOrder>::makeDerivative() const noexcept
requires (hasStaticCoefficients) {
    Derivative derivative;
    this->makeView().makeDerivative(derivative._coefficients);
    return derivative;
}


template <concepts::Numeric TValue, int PolynomialOrder>
std::span<const TValue> Polynomial<TValue,PolynomialOrder>::coefficients() const noexcept
requires (!hasStaticCoefficients) {
    return {_coefficients.data(), _coefficients.size()};
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr std::span<
    const TValue,
    Polynomial<TValue,PolynomialOrder>::coefficientCount
> Polynomial<TValue,PolynomialOrder>::coefficients() const noexcept
requires (hasStaticCoefficients) {
    return std::span<const TValue,coefficientCount>(
        _coefficients.data(),
        _coefficients.size());
}


template <concepts::Numeric TValue, int PolynomialOrder>
constexpr typename Polynomial<TValue,PolynomialOrder>::View
Polynomial<TValue,PolynomialOrder>::makeView() const noexcept {
    return View(_coefficients);
}


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TValue, int PolynomialOrder>
void GraphML::Serializer<maths::PolynomialView<TValue,PolynomialOrder>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class TValue, int PolynomialOrder>
void GraphML::Serializer<maths::PolynomialView<TValue,PolynomialOrder>>::operator()(Ref<GraphML::XMLElement> rElement,
                                                                                    Ref<const maths::PolynomialView<TValue,PolynomialOrder>> rInstance) {
    CIE_BEGIN_EXCEPTION_TRACING
    using SubSerializer = GraphML::Serializer<std::span<const TValue>>;
    SubSerializer subSerializer;
    GraphML::XMLElement child = rElement.addChild("polynomial");
    subSerializer(child, rInstance.coefficients());
    CIE_END_EXCEPTION_TRACING
}


template <class TValue, int PolynomialOrder>
void GraphML::Serializer<maths::Polynomial<TValue,PolynomialOrder>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class TValue, int PolynomialOrder>
void GraphML::Serializer<maths::Polynomial<TValue,PolynomialOrder>>::operator()(Ref<GraphML::XMLElement> rElement,
                                                                                Ref<const maths::Polynomial<TValue,PolynomialOrder>> rInstance) {
    CIE_BEGIN_EXCEPTION_TRACING
    using SubSerializer = GraphML::Serializer<std::span<const TValue>>;
    SubSerializer subSerializer;
    GraphML::XMLElement child = rElement.addChild("polynomial");
    subSerializer(child, rInstance.coefficients());
    CIE_END_EXCEPTION_TRACING
}


template <class TValue, int PolynomialOrder>
void GraphML::Deserializer<maths::Polynomial<TValue,PolynomialOrder>>::onElementBegin(Ptr<void> pThis,
                                                                                      std::string_view elementName,
                                                                                      [[maybe_unused]] std::span<GraphML::AttributePair> attributes) {
    CIE_BEGIN_EXCEPTION_TRACING

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    using SubDeserializer = GraphML::Deserializer<typename maths::Polynomial<TValue,PolynomialOrder>::Coefficients>;
    Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(rThis._coefficients, rThis.sax(), elementName);
    rThis.sax().push({
        pSubDeserializer,
        SubDeserializer::onElementBegin,
        SubDeserializer::onText,
        SubDeserializer::onElementEnd
    });

    //SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);

    CIE_END_EXCEPTION_TRACING
}


template <class TValue, int PolynomialOrder>
void GraphML::Deserializer<maths::Polynomial<TValue,PolynomialOrder>>::onText(Ptr<void>,
                                                                              std::string_view elementName) {
    CIE_THROW(
        Exception,
        "Unexpected text block while parsing a polynomial from <" << elementName << ">."
    )
}


template <class TValue, int PolynomialOrder>
void GraphML::Deserializer<maths::Polynomial<TValue,PolynomialOrder>>::onElementEnd(Ptr<void> pThis,
                                                                                    std::string_view elementName) {
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = maths::Polynomial<TValue,PolynomialOrder>(std::move(rThis._coefficients));
    rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem::io
