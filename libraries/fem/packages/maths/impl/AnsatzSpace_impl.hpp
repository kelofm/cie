#ifndef CIE_FEM_ANSATZ_SPACE_IMPL_HPP
#define CIE_FEM_ANSATZ_SPACE_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/maths/inc/power.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <class TScalarExpression, unsigned Dim>
inline void AnsatzSpaceDerivative<TScalarExpression,Dim>::evaluate(ConstIterator itArgumentBegin,
                                                                   ConstIterator itArgumentEnd,
                                                                   Iterator itOut) const
{
    const unsigned setSize = _ansatzSet.size();
    CIE_OUT_OF_RANGE_CHECK(setSize == _derivativeSet.size())
    CIE_OUT_OF_RANGE_CHECK(std::distance(itArgumentBegin, itArgumentEnd) == Dim)

    Ref<IndexBuffer> rIndexBuffer      = _buffer.template get<0>();
    Ref<ValueBuffer> rValueBuffer      = _buffer.template get<1>();
    Ref<ValueBuffer> rDerivativeBuffer = _buffer.template get<2>();

    // Fill the value and derivative buffers
    {
        Ptr<Value> pValue      = rValueBuffer.data();
        Ptr<Value> pDerivative = rDerivativeBuffer.data();

        for (auto it=itArgumentBegin; it!=itArgumentEnd; ++it) {
            const auto itEnd = it + 1;
            for (const auto& rScalarExpression : _ansatzSet) {
                rScalarExpression.evaluate(it, itEnd, pValue++);
            } // for scalarExpression in ansatzSet
            for (const auto& rScalarExpression : _derivativeSet) {
                rScalarExpression.evaluate(it, itEnd, pDerivative++);
            } // for scalarExpression in derivativeSet
        } // for component in arguments
    } // fill the value and derivative buffers

    // Compute the modified outer product
    for (unsigned iDerivative=0; iDerivative<Dim; ++iDerivative) {
        do {
            *itOut = static_cast<Value>(1);
            unsigned iIndex = 0;

            // First loop through the value buffer until
            // the derivative index is hit
            for (; iIndex<iDerivative; ++iIndex) {
                *itOut *= rValueBuffer[rIndexBuffer[iIndex] + iIndex * setSize];
            }

            // Then, use the derivative
            *itOut *= rDerivativeBuffer[rIndexBuffer[iIndex] + iIndex * setSize];

            // Finally, loop through the rest
            // of the value buffer
            for (++iIndex; iIndex<rIndexBuffer.size(); ++iIndex) {
                *itOut *= rValueBuffer[rIndexBuffer[iIndex] + iIndex * setSize];
            }

            ++itOut;
        } while (OuterProduct<Dim>::next(setSize, rIndexBuffer.begin()));
    } // for iDerivative in range(Dim)
}


template <class TScalarExpression, unsigned Dim>
unsigned AnsatzSpaceDerivative<TScalarExpression,Dim>::size() const noexcept
{
    return intPow(_ansatzSet.size(), Dim) * Dim;
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpaceDerivative<TScalarExpression,Dim>::AnsatzSpaceDerivative(Ref<const AnsatzSpace<TScalarExpression,Dim>> rAnsatzSpace)
    : _ansatzSet(rAnsatzSpace._set),
      _derivativeSet(),
      _buffer(IndexBuffer(),
              ValueBuffer(intPow(_ansatzSet.size(), Dim)),
              ValueBuffer(intPow(_ansatzSet.size(), Dim)))
{
    Ref<IndexBuffer> rIndexBuffer = _buffer.template get<0>();
    std::fill(rIndexBuffer.begin(),
              rIndexBuffer.end(),
              0u);

    _derivativeSet.reserve(_ansatzSet.size());
    for (const auto& rScalarExpression : _ansatzSet) {
        _derivativeSet.emplace_back(rScalarExpression.makeDerivative());
    }
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace() noexcept
    : AnsatzSpace(AnsatzSet {})
{
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace(AnsatzSet&& rSet) noexcept
    : _set(std::move(rSet)),
      _buffer(IndexBuffer(),
              ValueBuffer(intPow(_set.size(), Dim)))
{
    std::fill(_buffer.template get<0>().begin(),
              _buffer.template get<0>().end(),
              0u);
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace(const AnsatzSet& rSet)
    : AnsatzSpace(AnsatzSet(rSet))
{
}


template <class TScalarExpression, unsigned Dim>
inline void AnsatzSpace<TScalarExpression,Dim>::evaluate(ConstIterator itArgumentBegin,
                                                         ConstIterator itArgumentEnd,
                                                         Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(itArgumentBegin, itArgumentEnd) == Dim)
    Ref<IndexBuffer> rIndexBuffer = _buffer.template get<0>();
    Ref<ValueBuffer> rValueBuffer = _buffer.template get<1>();

    // No need to clear the index buffer of its leftover state
    // (the leftover state is all zeros)
    //std::fill(rIndexBuffer->begin(), rIndexBuffer->end(), 0);

    // Fill the value buffer
    Ptr<Value> pValue = rValueBuffer.data();
    for (; itArgumentBegin!=itArgumentEnd; itArgumentBegin++) {
        for (const auto& rScalarExpression : _set) {
            rScalarExpression.evaluate(itArgumentBegin,
                                        itArgumentBegin + 1,
                                        pValue++);
        } // for expression in expressions
    } // for argument in arguments

    const unsigned setSize = _set.size();
    do {
        // Compute product of bases
        *itOut = static_cast<Value>(1);
        for (unsigned iIndex=0; iIndex<rIndexBuffer.size(); ++iIndex) {
            *itOut *= rValueBuffer.at(rIndexBuffer.at(iIndex) + iIndex * setSize);
        }
        ++itOut;
    } while (OuterProduct<Dim>::next(setSize, rIndexBuffer.data()));
}


template <class TScalarExpression, unsigned Dim>
typename AnsatzSpace<TScalarExpression,Dim>::Derivative
AnsatzSpace<TScalarExpression,Dim>::makeDerivative() const
{
    return Derivative(*this);
}


template <class TScalarExpression, unsigned Dim>
unsigned AnsatzSpace<TScalarExpression,Dim>::size() const noexcept
{
    return intPow(_set.size(), Dim);
}


template <class TScalarExpression, unsigned Dim>
std::span<const TScalarExpression> AnsatzSpace<TScalarExpression,Dim>::ansatzSet() const noexcept
{
    return {_set.data(), _set.size()};
}


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TScalarExpression, unsigned Dim>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim>>::header(Ref<XMLElement> rElement)
{
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim>>::operator()(Ref<XMLElement> rElement,
                                                                                Ref<const maths::AnsatzSpace<TScalarExpression,Dim>> rInstance)
{
    CIE_BEGIN_EXCEPTION_TRACING

    static_assert(concepts::Container<std::span<const TScalarExpression>>, "test");
    using SubSerializer = GraphML::Serializer<std::span<const TScalarExpression>>;
    SubSerializer subSerializer;
    GraphML::XMLElement subElement = rElement.addChild("ansatz-space");
    subSerializer(subElement, rInstance.ansatzSet());
    //subSerializer(rElement, rInstance.ansatzSet());

    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim>>::onElementBegin(Ptr<void> pThis,
                                                                                      std::string_view elementName,
                                                                                      [[maybe_unused]] std::span<GraphML::AttributePair> attributes)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using SubDeserializer = GraphML::Deserializer<typename Value::AnsatzSet>;
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(rThis._set, rThis.sax(), elementName);

    rThis.sax().push({
        pSubDeserializer,
        SubDeserializer::onElementBegin,
        SubDeserializer::onText,
        SubDeserializer::onElementEnd
    });

    //SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);

    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim>>::onText(Ptr<void>,
                                                                              std::string_view)
{
    CIE_THROW(
        Exception,
        "Unexpected text block while parsing AnsatzSpace in GraphML."
    )
}


template <class TScalarExpression, unsigned Dim>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim>>::onElementEnd(Ptr<void> pThis,
                                                                                    std::string_view elementName)
{
    CIE_BEGIN_EXCEPTION_TRACING
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = Value(std::move(rThis._set));
    rThis.template release<Deserializer>(&rThis, elementName);
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem::io



#endif
