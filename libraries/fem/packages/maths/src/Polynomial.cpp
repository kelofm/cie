// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


template <concepts::Numeric TValue>
PolynomialView<TValue>::PolynomialView(ConstSpan coefficients) noexcept
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
Polynomial<TValue>& Polynomial<TValue>::operator=(const Polynomial& rRight) {
    _coefficients = rRight._coefficients;
    _wrapped = PolynomialView<TValue>(_coefficients);
    return *this;
}


template <class TValue>
Polynomial<TValue>::Polynomial(ConstSpan coefficients)
    : _coefficients(coefficients.begin(), coefficients.end())
{
    _wrapped = PolynomialView<TValue>(_coefficients);
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


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(PolynomialView);
CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(Polynomial);


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TValue>
void GraphML::Serializer<maths::PolynomialView<TValue>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
void GraphML::Serializer<maths::PolynomialView<TValue>>::operator()(Ref<GraphML::XMLElement> rElement,
                                                                    Ref<const maths::PolynomialView<TValue>> rInstance) {
    CIE_BEGIN_EXCEPTION_TRACING
    using SubSerializer = GraphML::Serializer<std::span<const TValue>>;
    SubSerializer subSerializer;
    GraphML::XMLElement child = rElement.addChild("polynomial");
    subSerializer(child, rInstance.coefficients());
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
void GraphML::Serializer<maths::Polynomial<TValue>>::header(Ref<GraphML::XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    defaultData.addAttribute("type", "polynomial");
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
void GraphML::Serializer<maths::Polynomial<TValue>>::operator()(Ref<GraphML::XMLElement> rElement,
                                                                Ref<const maths::Polynomial<TValue>> rInstance) {
    CIE_BEGIN_EXCEPTION_TRACING
    using SubSerializer = GraphML::Serializer<std::span<const TValue>>;
    SubSerializer subSerializer;
    GraphML::XMLElement child = rElement.addChild("polynomial");
    subSerializer(child, rInstance.coefficients());
    CIE_END_EXCEPTION_TRACING
}


template <class TValue>
void GraphML::Deserializer<maths::Polynomial<TValue>>::onElementBegin(Ptr<void> pThis,
                                                                      std::string_view elementName,
                                                                      [[maybe_unused]] std::span<GraphML::AttributePair> attributes) {
    CIE_BEGIN_EXCEPTION_TRACING

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    using SubDeserializer = GraphML::Deserializer<typename maths::Polynomial<TValue>::Coefficients>;
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


template <class TValue>
void GraphML::Deserializer<maths::Polynomial<TValue>>::onText(Ptr<void>,
                                                              std::string_view elementName) {
    CIE_THROW(
        Exception,
        "Unexpected text block while parsing a polynomial from <" << elementName << ">."
    )
}


template <class TValue>
void GraphML::Deserializer<maths::Polynomial<TValue>>::onElementEnd(Ptr<void> pThis,
                                                                    std::string_view elementName) {
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = maths::Polynomial<TValue>(std::move(rThis._coefficients));
    rThis.template release<Deserializer>(&rThis, elementName);
}


template struct GraphML::Serializer<maths::PolynomialView<float>>;
template struct GraphML::Serializer<maths::PolynomialView<double>>;
template struct GraphML::Serializer<maths::Polynomial<float>>;
template struct GraphML::Serializer<maths::Polynomial<double>>;
template struct GraphML::Deserializer<maths::Polynomial<float>>;
template struct GraphML::Deserializer<maths::Polynomial<double>>;


} // namespace cie::fem::io
