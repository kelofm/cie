#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/polynomial_evaluation.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem::maths {


template <concepts::Numeric T, int P>
constexpr PolynomialView<T,P>::PolynomialView() noexcept
    : _coefficients(static_cast<const T*>(nullptr), coefficientCount)
{}


template <concepts::Numeric T, int P>
PolynomialView<T,P>::PolynomialView(ConstSpan coefficients) noexcept
requires (!hasStaticCoefficients)
    : _coefficients(coefficients)
{}


template <concepts::Numeric T, int P>
constexpr PolynomialView<T,P>::PolynomialView(std::span<const T,coefficientCount> coefficients) noexcept
requires (hasStaticCoefficients)
    : _coefficients(coefficients)
{}


template <concepts::Numeric T, int P>
typename PolynomialView<T,P>::Derivative
PolynomialView<T,P>::makeDerivative(Span buffer) const
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
            T power = static_cast<T>(2);
            auto itBuffer = buffer.begin() + 1;
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power, ++itBuffer)
                *itBuffer = power * (*itCoefficient);
        } // if 2 < polynoialOrder
    } // if 1 < coefficientCount

    return PolynomialView(buffer);
}


template <concepts::Numeric T, int P>
constexpr typename PolynomialView<T,P>::Derivative
PolynomialView<T,P>::makeDerivative(std::span<T,Derivative::coefficientCount> buffer) const noexcept
requires (hasStaticCoefficients) {
    if constexpr (1 < coefficientCount) {
        // Push first coefficient (no multiplication required)
        buffer.front() = _coefficients[1];
        if constexpr (2 < coefficientCount) {
            const auto itCoefficientEnd = _coefficients.end();
            T power = static_cast<T>(2);
            auto itBuffer = buffer.begin() + 1;
            for (auto itCoefficient=_coefficients.begin()+2; itCoefficient!=itCoefficientEnd; ++itCoefficient, ++power, ++itBuffer)
                *itBuffer = power * (*itCoefficient);
        } // if 2 < polynoialOrder
    } /* if 1 < PolynomialOrder */

    return Derivative(buffer);
}


template <concepts::Numeric T, int P>
void PolynomialView<T,P>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan) const
requires (!hasStaticCoefficients) {
    out.front() = utils::evaluatePolynomial<utils::PolynomialEvaluation::Horner>(
        in.front(),
        _coefficients);
}


template <concepts::Numeric T, int P>
constexpr void PolynomialView<T,P>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan) const
requires (hasStaticCoefficients) {
    out.front() = utils::evaluatePolynomial<utils::PolynomialEvaluation::Horner>(
        in.front(),
        _coefficients);
}


template <concepts::Numeric T, int P>
constexpr unsigned PolynomialView<T,P>::size() noexcept {
    return 1u;
}


template <concepts::Numeric T, int P>
constexpr unsigned PolynomialView<T,P>::bufferSize() noexcept {
    return 0u;
}


template <concepts::Numeric T, int P>
typename PolynomialView<T,P>::ConstSpan
PolynomialView<T,P>::coefficients() const noexcept
requires (!hasStaticCoefficients) {
    return _coefficients;
}


template <concepts::Numeric T, int P>
constexpr std::span<const T,PolynomialView<T,P>::coefficientCount>
PolynomialView<T,P>::coefficients() const noexcept
requires (hasStaticCoefficients) {
    return _coefficients;
}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>::Polynomial(RightRef<Coefficients> rCoefficients) noexcept
    : _coefficients(std::move(rCoefficients))
{}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>::Polynomial(Polynomial&& rRight) noexcept
    : _coefficients(std::move(rRight._coefficients))
{}


template <concepts::Numeric T, int P, class TA>
Polynomial<T,P,TA>::Polynomial(Polynomial&& rRight)
requires (!hasStaticCoefficients)
    : _coefficients(std::move(rRight._coefficients))
{}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>::Polynomial(Polynomial&& rRight)
requires (hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric T, int P, class TA>
Polynomial<T,P,TA>::Polynomial(const Polynomial& rRight)
requires (!hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>::Polynomial(const Polynomial& rRight)
requires (hasStaticCoefficients)
    : _coefficients(rRight._coefficients)
{}


template <concepts::Numeric T, int P, class TA>
Polynomial<T,P,TA>::Polynomial(ConstSpan coefficients)
requires (!hasStaticCoefficients)
    : _coefficients(coefficients.begin(), coefficients.end())
{}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>::Polynomial(std::span<const T,coefficientCount> coefficients)
requires (hasStaticCoefficients)
    : _coefficients() {
        std::copy_n(
            coefficients.data(),
            coefficientCount,
            _coefficients.data());
}


template <concepts::Numeric T, int P, class TA>
Polynomial<T,P,TA>&
Polynomial<T,P,TA>::operator=(Polynomial&& rRight) noexcept
requires (!hasStaticCoefficients) {
    _coefficients = std::move(rRight._coefficients);
    return *this;
}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>&
Polynomial<T,P,TA>::operator=(Polynomial&& rRight) noexcept
requires (hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric T, int P, class TA>
Polynomial<T,P,TA>& Polynomial<T,P,TA>::operator=(const Polynomial& rRight)
requires (!hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric T, int P, class TA>
constexpr Polynomial<T,P,TA>& Polynomial<T,P,TA>::operator=(const Polynomial& rRight) noexcept
requires (hasStaticCoefficients) {
    _coefficients = rRight._coefficients;
    return *this;
}


template <concepts::Numeric T, int P, class TA>
void Polynomial<T,P,TA>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const
requires (!hasStaticCoefficients) {
    this->makeView().evaluate(in, out, buffer);
}


template <concepts::Numeric T, int P, class TA>
constexpr void Polynomial<T,P,TA>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const noexcept
requires (hasStaticCoefficients) {
    this->makeView().evaluate(in, out, buffer);
}


template <concepts::Numeric T, int P, class TA>
constexpr unsigned Polynomial<T,P,TA>::size() noexcept {
    return View::size();
}


template <concepts::Numeric T, int P, class TA>
constexpr unsigned Polynomial<T,P,TA>::bufferSize() noexcept {
    return View::bufferSize();
}


template <concepts::Numeric T, int P, class TA>
typename Polynomial<T,P,TA>::Derivative
Polynomial<T,P,TA>::makeDerivative() const
requires (!hasStaticCoefficients) {
    Polynomial derivative;
    derivative._coefficients.resize(_coefficients.empty() ? 0 : _coefficients.size() - 1);
    this->makeView().makeDerivative(derivative._coefficients);
    return derivative;
}


template <concepts::Numeric T, int P, class TA>
constexpr typename Polynomial<T,P,TA>::Derivative
Polynomial<T,P,TA>::makeDerivative() const noexcept
requires (hasStaticCoefficients) {
    Derivative derivative;
    this->makeView().makeDerivative(derivative._coefficients);
    return derivative;
}


template <concepts::Numeric T, int P, class TA>
std::span<const T> Polynomial<T,P,TA>::coefficients() const noexcept
requires (!hasStaticCoefficients) {
    return {_coefficients.data(), _coefficients.size()};
}


template <concepts::Numeric T, int P, class TA>
constexpr std::span<
    const T,
    Polynomial<T,P,TA>::coefficientCount
> Polynomial<T,P,TA>::coefficients() const noexcept
requires (hasStaticCoefficients) {
    return std::span<const T,coefficientCount>(
        _coefficients.data(),
        _coefficients.size());
}


template <concepts::Numeric T, int P, class TA>
constexpr typename Polynomial<T,P,TA>::View
Polynomial<T,P,TA>::makeView() const noexcept {
    return View(_coefficients);
}


template <concepts::Numeric T, int P, class TA>
void Polynomial<T,P,TA>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        if constexpr (!hasStaticCoefficients)
            cie::io::BinarySerializer::serialize(
                rStream,
                static_cast<std::size_t>(_coefficients.size()));
        cie::io::BinarySerializer::serialize(
            rStream,
            _coefficients.data(),
            _coefficients.size());
}


template <concepts::Numeric T, int P, class TA>
void Polynomial<T,P,TA>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<Polynomial> rInstance,
    tags::Binary) {
        std::size_t coefficientCount = static_cast<std::size_t>(P);
        if constexpr (!hasStaticCoefficients) {
            cie::io::BinarySerializer::deserialize(
                rStream,
                coefficientCount);
            rInstance._coefficients.resize(coefficientCount);
        }
        cie::io::BinarySerializer::deserialize(
            rStream,
            rInstance._coefficients.data(),
            rInstance._coefficients.size());
}


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class T, int P>
void GraphML::Serializer<maths::PolynomialView<T,P>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
        GraphML::XMLElement defaultData = rElement.addChild("default");
        defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class T, int P>
void GraphML::Serializer<maths::PolynomialView<T,P>>::operator()(
    Ref<GraphML::XMLElement> rElement,
    Ref<const maths::PolynomialView<T,P>> rInstance) {
        CIE_BEGIN_EXCEPTION_TRACING
            using SubSerializer = GraphML::Serializer<std::span<const T>>;
            SubSerializer subSerializer;
            GraphML::XMLElement child = rElement.addChild("polynomial");
            subSerializer(child, rInstance.coefficients());
        CIE_END_EXCEPTION_TRACING
}


template <class T, int P, class TA>
void GraphML::Serializer<maths::Polynomial<T,P,TA>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
        GraphML::XMLElement defaultData = rElement.addChild("default");
        defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class T, int P, class TA>
void GraphML::Serializer<maths::Polynomial<T,P,TA>>::operator()(
    Ref<GraphML::XMLElement> rElement,
    Ref<const maths::Polynomial<T,P,TA>> rInstance) {
        CIE_BEGIN_EXCEPTION_TRACING
            using SubSerializer = GraphML::Serializer<std::span<const T>>;
            SubSerializer subSerializer;
            GraphML::XMLElement child = rElement.addChild("polynomial");
            subSerializer(child, rInstance.coefficients());
        CIE_END_EXCEPTION_TRACING
}


template <class T, int P, class TA>
void GraphML::Deserializer<maths::Polynomial<T,P,TA>>::onElementBegin(
    Ptr<void> pThis,
    std::string_view elementName,
    [[maybe_unused]] std::span<GraphML::AttributePair> attributes) {
        CIE_BEGIN_EXCEPTION_TRACING
            Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
            using SubDeserializer = GraphML::Deserializer<typename maths::Polynomial<T,P,TA>::Coefficients>;
            Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(rThis._coefficients, rThis.sax(), elementName);
            rThis.sax().push({
                pSubDeserializer,
                SubDeserializer::onElementBegin,
                SubDeserializer::onText,
                SubDeserializer::onElementEnd});
        CIE_END_EXCEPTION_TRACING
}


template <class T, int P, class TA>
void GraphML::Deserializer<maths::Polynomial<T,P,TA>>::onText(
    Ptr<void>,
    std::string_view elementName) {
        CIE_THROW(
            Exception,
            "Unexpected text block while parsing a polynomial from <" << elementName << ">."
        )
}


template <class T, int P, class TA>
void GraphML::Deserializer<maths::Polynomial<T,P,TA>>::onElementEnd(
    Ptr<void> pThis,
    std::string_view elementName) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.instance() = maths::Polynomial<T,P,TA>(std::move(rThis._coefficients));
        rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem::io
