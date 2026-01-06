#pragma once

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
BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::BufferedAnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet,
                                                                                    std::span<const typename TScalarExpression::Derivative> derivativeSet)
    : _ansatzSet(ansatzSet),
      _derivativeSet(derivativeSet),
      _buffer()
{
    if (_derivativeSet.size() != _ansatzSet.size()) {
        CIE_THROW(
            OutOfRangeException,
               "number of provided ansatz derivatives (" << _derivativeSet.size() << ") "
            << "does not match the number of ansatz functions (" << _ansatzSet.size() << ")"
        )
    }
}


template <class TScalarExpression, unsigned Dim>
BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::BufferedAnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet,
                                                                                    std::span<const typename TScalarExpression::Derivative> derivativeSet,
                                                                                    Span buffer)
    : BufferedAnsatzSpaceDerivative(ansatzSet, derivativeSet)
{
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim>
void BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::evaluate(ConstSpan in, Span out) const {
    const unsigned setSize = _ansatzSet.size();
    CIE_OUT_OF_RANGE_CHECK(in.size() == Dim)
    CIE_OUT_OF_RANGE_CHECK(setSize == _derivativeSet.size())
    CIE_OUT_OF_RANGE_CHECK(out.size() == this->size())

    auto indexBuffer        = this->getIndexBuffer();
    auto ansatzBuffer       = this->getAnsatzBuffer();
    auto derivativeBuffer   = this->getDerivativeBuffer();

    // Fill the value and derivative buffers
    {
        Ptr<Value> pValue      = ansatzBuffer.data();
        Ptr<Value> pDerivative = derivativeBuffer.data();

        for (auto c : in) {
            ConstSpan scalarIn {&c, 1};
            for (auto& rScalarExpression : _ansatzSet) {
                rScalarExpression.evaluate(scalarIn, {pValue, pValue + 1});
                ++pValue;
            } // for scalarExpression in ansatzSet
            for (auto& rScalarExpression : _derivativeSet) {
                rScalarExpression.evaluate(scalarIn, {pDerivative, pDerivative + 1});
                ++pDerivative;
            } // for scalarExpression in derivativeSet
        } // for component in arguments
    } // fill the value and derivative buffers

    // Compute the modified outer product
    auto itOut = out.data();
    for (unsigned iDerivative=0; iDerivative<Dim; ++iDerivative) {
        do {
            *itOut = static_cast<Value>(1);
            unsigned iIndex = 0;

            // First loop through the value buffer until
            // the derivative index is hit
            for (; iIndex<iDerivative; ++iIndex) {
                *itOut *= ansatzBuffer[indexBuffer[iIndex] + iIndex * setSize];
            }

            // Then, use the derivative
            *itOut *= derivativeBuffer[indexBuffer[iIndex] + iIndex * setSize];

            // Finally, loop through the rest
            // of the value buffer
            for (++iIndex; iIndex<indexBuffer.size(); ++iIndex) {
                *itOut *= ansatzBuffer[indexBuffer[iIndex] + iIndex * setSize];
            }

            ++itOut;
        } while (cie::maths::OuterProduct<Dim>::next(setSize, indexBuffer.data()));
    } // for iDerivative in range(Dim)
}


template <class TScalarExpression, unsigned Dim>
unsigned BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::size() const noexcept {
    return intPow(_ansatzSet.size(), Dim) * Dim;
}


template <class TScalarExpression, unsigned Dim>
unsigned BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::getMinBufferSize() const noexcept {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the ansatz buffer from the begin of the buffer, in bytes.
    constexpr unsigned ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    const unsigned valueBufferSize = intPow(_ansatzSet.size(), Dim) * sizeof(Value);
    const unsigned derivativeBufferOffset = ansatzBufferOffset + valueBufferSize;

    return (derivativeBufferOffset + valueBufferSize) / sizeof(Value);
}


template <class TScalarExpression, unsigned Dim>
void BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::setBuffer(Span buffer) {
    if (buffer.size() < this->getMinBufferSize()) {
        CIE_THROW(
            OutOfRangeException ,
               "provided buffer size (" << buffer.size() << ") "
            << "does not meet the minimum requirement (" << this->getMinBufferSize() << ")"
        )
    }

    _buffer = buffer;

    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim>
std::span<unsigned,Dim> BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::getIndexBuffer() const noexcept {
    return std::span<unsigned,Dim> {
        static_cast<unsigned*>(static_cast<void*>(_buffer.data())),
        Dim
    };
}


template <class TScalarExpression, unsigned Dim>
typename BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::Span
BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::getAnsatzBuffer() const noexcept {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    return Span(
        _buffer.data() + ansatzBufferOffset / sizeof(Value),
        intPow(_ansatzSet.size(), Dim)
    );
}


template <class TScalarExpression, unsigned Dim>
typename BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::Span
BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::getDerivativeBuffer() const noexcept {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned valueBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);
    const unsigned valueBufferSize = intPow(_ansatzSet.size(), Dim) * sizeof(Value);

    return Span(
        _buffer.data() + (valueBufferOffset + valueBufferSize) / sizeof(Value),
        intPow(_ansatzSet.size(), Dim)
    );
}


template <class TScalarExpression, unsigned Dim>
std::span<const TScalarExpression>
BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>::ansatzSet() const noexcept {
    return _ansatzSet;
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpaceDerivative<TScalarExpression,Dim>::AnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet)
    : _ansatzSet(ansatzSet.begin(), ansatzSet.end()),
      _derivativeSet(ansatzSet.size()),
      _buffer(),
      _wrapped()
{
    std::transform(
        _ansatzSet.begin(),
        _ansatzSet.end(),
        _derivativeSet.begin(),
        [](Ref<const TScalarExpression> rAnsatzFunction){
             return rAnsatzFunction.makeDerivative();
        });
    _wrapped = BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>(
        {_ansatzSet.data(), _ansatzSet.size()},
        {_derivativeSet.data(), _derivativeSet.size()}
    );
    _buffer.resize(_wrapped.getMinBufferSize());
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpaceDerivative<TScalarExpression,Dim>::AnsatzSpaceDerivative(const AnsatzSpaceDerivative& rRhs)
    : AnsatzSpaceDerivative({rRhs._set.data(), rRhs._set.size()})
{
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpaceDerivative<TScalarExpression,Dim>&
AnsatzSpaceDerivative<TScalarExpression,Dim>::operator=(const AnsatzSpaceDerivative& rRhs) {
    (*this) = AnsatzSpaceDerivative(rRhs);
}


template <class TScalarExpression, unsigned Dim>
void AnsatzSpaceDerivative<TScalarExpression,Dim>::evaluate(ConstSpan in, Span out) const {
    _wrapped.evaluate(in, out);
}


template <class TScalarExpression, unsigned Dim>
unsigned AnsatzSpaceDerivative<TScalarExpression,Dim>::size() const noexcept {
    return _wrapped.size();
}


template <class TScalarExpression, unsigned Dim>
typename AnsatzSpaceDerivative<TScalarExpression,Dim>::Buffered
AnsatzSpaceDerivative<TScalarExpression,Dim>::makeBuffered() const noexcept {
    return BufferedAnsatzSpaceDerivative<TScalarExpression,Dim>(
        _ansatzSet,
        _derivativeSet);
}


template <class TScalarExpression, unsigned Dim>
BufferedAnsatzSpace<TScalarExpression,Dim>::BufferedAnsatzSpace() noexcept
    : _set(),
      _buffer()
{
}


template <class TScalarExpression, unsigned Dim>
BufferedAnsatzSpace<TScalarExpression,Dim>::BufferedAnsatzSpace(std::span<const TScalarExpression> ansatzSet) noexcept
    : _set(ansatzSet),
      _buffer()
{
}


template <class TScalarExpression, unsigned Dim>
BufferedAnsatzSpace<TScalarExpression,Dim>::BufferedAnsatzSpace(std::span<const TScalarExpression> ansatzSet,
                                                                Span buffer)
    : _set(ansatzSet),
      _buffer() {
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim>
void BufferedAnsatzSpace<TScalarExpression,Dim>::evaluate(ConstSpan in, Span out) const {
    // Sanity checks
    CIE_OUT_OF_RANGE_CHECK(in.size() == Dim)
    CIE_OUT_OF_RANGE_CHECK(out.size() == this->size())

    auto valueBuffer = this->getValueBuffer();
    auto indexBuffer = this->getIndexBuffer();

    // No need to clear the index buffer of its leftover state
    // (the leftover state is all zeros)
    //std::fill(indexBuffer->begin(), rIndexBuffer->end(), 0);

    // Fill the value buffer
    auto pValue = valueBuffer.data();
    for (unsigned iDim=0u; iDim<Dim; ++iDim) {
        for (auto& rScalarExpression : _set) {
            rScalarExpression.evaluate({in.data() + iDim, 1}, {pValue++, 1});
        } // for rScalarExpression in _set
    } // for iDim in range(Dim)

    const unsigned setSize = _set.size();
    auto itOut = out.data();
    do {
        // Compute product of bases
        *itOut = static_cast<Value>(1);
        for (unsigned iIndex=0; iIndex<indexBuffer.size(); ++iIndex) {
            *itOut *= valueBuffer[indexBuffer[iIndex] + iIndex * setSize];
        }
        ++itOut;
    } while (cie::maths::OuterProduct<Dim>::next(setSize, indexBuffer.data()));
}


template <class TScalarExpression, unsigned Dim>
unsigned BufferedAnsatzSpace<TScalarExpression,Dim>::size() const noexcept {
    return intPow(_set.size(), Dim);
}


template <class TScalarExpression, unsigned Dim>
unsigned BufferedAnsatzSpace<TScalarExpression,Dim>::getMinBufferSize() const noexcept {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned valueBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    // Required size of the value buffer in bytes.
    unsigned valueBufferSize = sizeof(Value) * intPow(_set.size(), Dim);

    return (valueBufferOffset + valueBufferSize) / sizeof(Value);
}


template <class TScalarExpression, unsigned Dim>
void BufferedAnsatzSpace<TScalarExpression,Dim>::setBuffer(Span buffer) {
    if (buffer.size() < this->getMinBufferSize()) {
        CIE_THROW(
            OutOfRangeException ,
               "provided buffer size (" << buffer.size() << ") "
            << "does not meet the minimum requirement (" << this->getMinBufferSize() << ")"
        )
    }

    _buffer = buffer;

    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim>
std::span<unsigned,Dim> BufferedAnsatzSpace<TScalarExpression,Dim>::getIndexBuffer() const noexcept {
    return std::span<unsigned,Dim> {
        static_cast<unsigned*>(static_cast<void*>(_buffer.data())),
        Dim
    };
}


template <class TScalarExpression, unsigned Dim>
typename BufferedAnsatzSpace<TScalarExpression,Dim>::Span
BufferedAnsatzSpace<TScalarExpression,Dim>::getValueBuffer() const noexcept {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned valueBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    return Span(
        _buffer.data() + valueBufferOffset / sizeof(Value),
        intPow(_set.size(), Dim)
    );
}


template <class TScalarExpression, unsigned Dim>
std::span<const TScalarExpression> BufferedAnsatzSpace<TScalarExpression,Dim>::ansatzSet() const noexcept {
    return _set;
}


template <class TScalarExpression, unsigned Dim>
typename AnsatzSpace<TScalarExpression,Dim>::Buffered
AnsatzSpace<TScalarExpression,Dim>::makeBuffered() const noexcept {
    return BufferedAnsatzSpace<TScalarExpression,Dim>(_set);
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace() noexcept
    : AnsatzSpace(AnsatzSet {})
{
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace(AnsatzSet&& rSet) noexcept
    : _set(std::move(rSet)),
      _buffer(),
      _wrapped()
{
    _wrapped = BufferedAnsatzSpace<TScalarExpression,Dim>({_set.data(), _set.size()});
    _buffer.resize(_wrapped.getMinBufferSize());
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace(const AnsatzSet& rSet)
    : AnsatzSpace(AnsatzSet(rSet))
{
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>::AnsatzSpace(const AnsatzSpace& rRhs)
    : AnsatzSpace(rRhs._set)
{
}


template <class TScalarExpression, unsigned Dim>
AnsatzSpace<TScalarExpression,Dim>&
AnsatzSpace<TScalarExpression,Dim>::operator=(const AnsatzSpace& rRhs) {
    (*this) = AnsatzSpace(rRhs._set);
    return *this;
}


template <class TScalarExpression, unsigned Dim>
void AnsatzSpace<TScalarExpression,Dim>::evaluate(ConstSpan in, Span out) const {
    _wrapped.evaluate(in, out);
}


template <class TScalarExpression, unsigned Dim>
typename AnsatzSpace<TScalarExpression,Dim>::Derivative
AnsatzSpace<TScalarExpression,Dim>::makeDerivative() const {
    return Derivative(this->ansatzSet());
}


template <class TScalarExpression, unsigned Dim>
unsigned AnsatzSpace<TScalarExpression,Dim>::size() const noexcept {
    return _wrapped.size();
}


template <class TScalarExpression, unsigned Dim>
std::span<const TScalarExpression> AnsatzSpace<TScalarExpression,Dim>::ansatzSet() const noexcept {
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
