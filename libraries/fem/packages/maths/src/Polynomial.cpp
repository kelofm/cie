// --- FEM Includes ---
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem::maths {


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
