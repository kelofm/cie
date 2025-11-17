#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/types/inc/tags.hpp" // tags::Binary

// --- STL Includes ---
#include <algorithm> // std::copy, std::transform (already included anyway)
#include <functional> // std::multiplies
#include <numeric> // std::accumulate


namespace cie::fem::maths {


template <concepts::Numeric TValue, unsigned Dimension>
void ScaleTranslateTransformDerivative<TValue,Dimension>::evaluate(ConstSpan, Span out) const noexcept
{
    // Return a Dimension x Dimension matrix with _scales on the main diagonal.
    static_assert(0 < Dimension);
    auto itScale = this->_scales.cbegin();
    auto itOut = out.data();

    *itOut++ = *itScale++;
    for (unsigned iColumn=0; iColumn<Dimension-1; ++iColumn) {
        for (unsigned iNullComponent=0; iNullComponent<Dimension; ++iNullComponent) {
            *itOut++ = static_cast<TValue>(0);
        } // for iNullComponent in range(Dimension)
        *itOut++ = *itScale++;
    } // for iColumn in range(Dimension-1)
}


template <concepts::Numeric TValue, unsigned Dimension>
TValue ScaleTranslateTransformDerivative<TValue,Dimension>::evaluateDeterminant(ConstSpan) const noexcept
{
    return std::accumulate(this->_scales.begin(),
                           this->_scales.end(),
                           static_cast<TValue>(1),
                           std::multiplies<TValue>());
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned ScaleTranslateTransformDerivative<TValue,Dimension>::size() const noexcept
{
    return Dimension * Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                             Ref<ScaleTranslateTransformDerivative<TValue,Dimension>> rObject)
{
    for (const auto scale : rObject._scales)
        rStream << scale << ' ';
    return rStream;
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform(TPointIt itTransformedBegin,
                                                                   [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(
        std::distance(itTransformedBegin, itTransformedEnd) == 2,
        "Expecting 2 points, but got " << std::distance(itTransformedBegin, itTransformedEnd)
    )

    const auto& rBase = *itTransformedBegin;
    const auto& rOp = *(itTransformedBegin + 1);
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        const TValue diff = rOp[iDim] - rBase[iDim];
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[iDim] = diff / static_cast<TValue>(2);
        this->_offset[iDim] = (rOp[iDim] + rBase[iDim]) / static_cast<TValue>(2);
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
void ScaleTranslateTransform<TValue,Dimension>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(in.size() == Dimension)
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        out[iDim] = in[iDim] * this->_scales[iDim] + this->_offset[iDim];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
ScaleTranslateTransform<TValue,Dimension>::ScaleTranslateTransform() noexcept
{
    std::fill(this->_scales.begin(),
              this->_scales.end(),
              static_cast<TValue>(1));
}


template <concepts::Numeric TValue, unsigned Dimension>
template <concepts::Iterator TPointIt>
TranslateScaleTransform<TValue,Dimension>::TranslateScaleTransform(TPointIt itTransformedBegin,
                                                                   [[maybe_unused]] TPointIt itTransformedEnd)
{
    CIE_OUT_OF_RANGE_CHECK(
        std::distance(itTransformedBegin, itTransformedEnd) == 2,
        "Expecting 2 points, but got " << std::distance(itTransformedBegin, itTransformedEnd)
    )

    const auto& rBase = *itTransformedBegin;
    const auto& rOp = *(itTransformedBegin + 1);
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        const TValue diff = (rOp[iDim] - rBase[iDim]);
        CIE_DIVISION_BY_ZERO_CHECK(std::numeric_limits<TValue>::epsilon() < std::abs(diff))
        this->_scales[iDim] = diff / static_cast<TValue>(2);
        this->_offset[iDim] = (rOp[iDim] + rBase[iDim]) / diff;
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
void TranslateScaleTransform<TValue,Dimension>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(in.size() == Dimension)
    for (unsigned iDim=0; iDim<Dimension; ++iDim) {
        out[iDim] = (in[iDim] + this->_offset[iDim]) * this->_scales[iDim];
    }
}


template <concepts::Numeric TValue, unsigned Dimension>
unsigned TranslateScaleTransform<TValue,Dimension>::size() const noexcept
{
    return Dimension;
}


template <concepts::Numeric TValue, unsigned Dimension>
Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                             Ref<const ScaleTranslateTransform<TValue,Dimension>> rObject)
{
    StaticArray<TValue,Dimension> input, output;

    std::fill(input.begin(), input.end(), static_cast<TValue>(-1));
    rObject.evaluate(input, output);
    for (const auto c : output) rStream << c << ' ';

    std::fill(input.begin(), input.end(), static_cast<TValue>(1));
    rObject.evaluate(input, output);
    for (const auto c : output) rStream << c << ' ';

    return rStream;
}


template <concepts::Numeric TValue, unsigned Dimension>
Ref<std::ostream> operator<<(Ref<std::ostream> rStream,
                             Ref<const TranslateScaleTransform<TValue,Dimension>> rObject)
{
    for (const auto scale : rObject._scales) rStream << scale << ' ';
    for (const auto component : rObject._offset) rStream << component << ' ';
    return rStream;
}


} // namespace cie::fem::maths


// --- IO --- //


namespace cie::fem::io {


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::ScaleTranslateTransformDerivative<TValue,Dimension>>::header(Ref<XMLElement> rElement) noexcept
{
    auto descriptionElement = rElement.addChild("desc");
    std::stringstream description;
    description << "Derivative of a scaling in " << Dimension << " dimensions.";
    descriptionElement.setValue(description.view());
}


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::ScaleTranslateTransformDerivative<TValue,Dimension>>::operator()(Ref<XMLElement> rElement,
                                                                                                 Ref<const maths::ScaleTranslateTransformDerivative<TValue,Dimension>> rObject) noexcept
{
    std::stringstream stream;
    stream << rObject;
    GraphML::XMLElement subElement = rElement.addChild("d-st-tr");
    subElement.setValue(stream.view());
}


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::ScaleTranslateTransform<TValue,Dimension>>::header(Ref<XMLElement> rElement) noexcept
{
    auto descriptionElement = rElement.addChild("desc");
    std::stringstream description;
    description << "Scaling followed by a translation in " << Dimension << " dimensions.";
    descriptionElement.setValue(description.view());
}


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::ScaleTranslateTransform<TValue,Dimension>>::operator()(Ref<XMLElement> rElement,
                                                                                       Ref<const maths::ScaleTranslateTransform<TValue,Dimension>> rObject) noexcept
{
    StaticArray<TValue,Dimension>   input;
    StaticArray<TValue,2*Dimension> output;

    std::fill(input.begin(), input.end(), static_cast<TValue>(-1));
    rObject.evaluate(input, {output.data(), output.data() + Dimension});

    std::fill(input.begin(), input.end(), static_cast<TValue>(1));
    rObject.evaluate(input, {output.data() + Dimension, output.data() + 2 * Dimension});

    GraphML::XMLElement child = rElement.addChild("st-tr");
    using SubSerializer = GraphML::Serializer<std::span<const TValue>>;
    SubSerializer subSerializer;
    //subSerializer.setFormat(tags::Binary::flags());

    CIE_BEGIN_EXCEPTION_TRACING
    subSerializer(child, std::span<const TValue>(output));
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric TValue, unsigned Dimension>
void io::GraphML::Deserializer<maths::ScaleTranslateTransform<TValue,Dimension>>::onElementBegin(Ptr<void> pThis,
                                                                                                 std::string_view elementName,
                                                                                                 std::span<GraphML::AttributePair>)
{
    if (elementName != "st-tr") {
        CIE_THROW(
            Exception,
            "Expecting an <st-tr> element but got <" << elementName << ">."
        )
    }

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    using SubDeserializer = io::GraphML::Deserializer<decltype(Deserializer::_buffer)>;
    Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(rThis._buffer, rThis.sax(), elementName);
    //pSubDeserializer->setFormat(tags::Binary::flags());

    rThis.sax().push({
        pSubDeserializer,
        SubDeserializer::onElementBegin,
        SubDeserializer::onText,
        SubDeserializer::onElementEnd
    });
}


template <concepts::Numeric TValue, unsigned Dimension>
void io::GraphML::Deserializer<maths::ScaleTranslateTransform<TValue,Dimension>>::onText(Ptr<void>,
                                                                                         std::string_view)
{
    CIE_THROW(
        Exception,
        "Unexpected text block while parsing <st-tr>."
    )
}


template <concepts::Numeric TValue, unsigned Dimension>
void io::GraphML::Deserializer<maths::ScaleTranslateTransform<TValue,Dimension>>::onElementEnd(Ptr<void> pThis,
                                                                                               std::string_view elementName)
{
    if (elementName != "st-tr") {
        CIE_THROW(
            Exception,
            "Expecting a closing tag for <st-tr> but got one for <" << elementName << ">."
        )
    }
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    StaticArray<std::span<TValue>,2> transformed {
        std::span<TValue>(rThis._buffer.data(),             rThis._buffer.data() + Dimension),
        std::span<TValue>(rThis._buffer.data() + Dimension, rThis._buffer.data() + 2 * Dimension)
    };
    rThis.instance() = maths::ScaleTranslateTransform<TValue,Dimension>(transformed.begin(), transformed.end());

    rThis.template release<Deserializer>(&rThis, elementName);
}


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::TranslateScaleTransform<TValue,Dimension>>::header(Ref<XMLElement> rElement) noexcept
{
    auto descriptionElement = rElement.addChild("desc");
    std::stringstream description;
    description << "Translation followed by scaling in " << Dimension << " dimensions.";
    descriptionElement.setValue(description.view());
}


template <concepts::Numeric TValue, unsigned Dimension>
void GraphML::Serializer<maths::TranslateScaleTransform<TValue,Dimension>>::operator()(Ref<XMLElement> rElement,
                                                                                       Ref<const maths::TranslateScaleTransform<TValue,Dimension>> rObject) noexcept
{
    std::stringstream stream;
    stream << rObject;
    GraphML::XMLElement subElement = rElement.addChild("ts-tr");
    subElement.setValue(stream.view());
}


} // namespace cie::fem::io
