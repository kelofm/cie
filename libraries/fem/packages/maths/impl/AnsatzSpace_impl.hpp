#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression,SetSize> ansatzSet,
    std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet) noexcept
requires hasStaticBasis
    : _ansatzSet(ansatzSet),
      _derivativeSet(derivativeSet),
      _buffer()
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression,SetSize> ansatzSet,
    std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet,
    std::span<Value,staticBufferSize> buffer) noexcept
requires hasStaticBasis
    : AnsatzSpaceDerivativeView(ansatzSet, derivativeSet)
{
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression> ansatzSet,
    std::span<const typename TScalarExpression::Derivative> derivativeSet)
requires (!hasStaticBasis)
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


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression> ansatzSet,
    std::span<const typename TScalarExpression::Derivative> derivativeSet,
    Span buffer)
requires (!hasStaticBasis)
    : AnsatzSpaceDerivativeView(ansatzSet, derivativeSet)
{
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::evaluate(ConstSpan in, Span out) const {
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
                rScalarExpression.evaluate(scalarIn, {pValue, 1});
                ++pValue;
            } // for scalarExpression in ansatzSet
            for (auto& rScalarExpression : _derivativeSet) {
                rScalarExpression.evaluate(scalarIn, {pDerivative, 1});
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


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return intPow(_ansatzSet.size(), Dim) * Dim;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::size() noexcept
requires hasStaticBasis {
    return intPow(SetSize, Dim) * Dim;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getMinBufferSize() const noexcept
requires (!hasStaticBasis) {
    const unsigned valueBufferSize = intPow(_ansatzSet.size(), Dim) * sizeof(Value);
    const unsigned derivativeBufferOffset = ansatzBufferOffset + valueBufferSize;
    return (derivativeBufferOffset + valueBufferSize) / sizeof(Value);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getMinBufferSize() noexcept
requires hasStaticBasis {
    return staticBufferSize;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::setBuffer(Span buffer)
requires (!hasStaticBasis) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr void AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::setBuffer(std::span<Value,staticBufferSize> buffer) noexcept
requires hasStaticBasis {
    _buffer = buffer;
    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<unsigned,Dim> AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getIndexBuffer() const noexcept {
    return std::span<unsigned,Dim> {
        static_cast<unsigned*>(static_cast<void*>(_buffer.data())),
        Dim
    };
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::Span
AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getAnsatzBuffer() const noexcept
requires (!hasStaticBasis) {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    return Span(
        _buffer.data() + ansatzBufferOffset / sizeof(Value),
        intPow(_ansatzSet.size(), Dim)
    );
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<
    typename AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::Value,
    AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::valueBufferSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getAnsatzBuffer() const noexcept
requires hasStaticBasis {
    return Span(
        _buffer.data() + ansatzBufferOffset / sizeof(Value),
        intPow(SetSize, Dim)
    );
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::Span
AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getDerivativeBuffer() const noexcept
requires (!hasStaticBasis) {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);
    const unsigned valueBufferSize = intPow(_ansatzSet.size(), Dim) * sizeof(Value);

    return Span(
        _buffer.data() + (ansatzBufferOffset + valueBufferSize) / sizeof(Value),
        intPow(_ansatzSet.size(), Dim)
    );
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<
    typename AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::Value,
    AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::valueBufferSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::getDerivativeBuffer() const noexcept
requires hasStaticBasis {
    return Span(
        _buffer.data() + (ansatzBufferOffset + valueBufferSize) / sizeof(Value),
        intPow(SetSize, Dim)
    );
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept {
    return _ansatzSet;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const typename TScalarExpression::Derivative>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::derivativeSet() const noexcept {
    return _derivativeSet;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(
    std::span<const TScalarExpression,SetSize> ansatzSet) noexcept
requires hasStaticBasis {
    std::copy_n(
        ansatzSet.data(),
        SetSize,
        _ansatzSet.data());
    std::transform(
        _ansatzSet.begin(),
        _ansatzSet.end(),
        _derivativeSet.begin(),
        [](Ref<const TScalarExpression> rAnsatzFunction){
             return rAnsatzFunction.makeDerivative();
        });
    _wrapped = AnsatzSpaceDerivativeView<TScalarExpression,Dim>(
        {_ansatzSet.data(), _ansatzSet.size()},
        {_derivativeSet.data(), _derivativeSet.size()});
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet)
requires (!hasStaticBasis)
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
    _wrapped = AnsatzSpaceDerivativeView<TScalarExpression,Dim>(
        {_ansatzSet.data(), _ansatzSet.size()},
        {_derivativeSet.data(), _derivativeSet.size()});
    _buffer.resize(_wrapped.getMinBufferSize());
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(const AnsatzSpaceDerivative& rRhs)
requires (!hasStaticBasis)
    : _ansatzSet(rRhs._ansatzSet),
      _derivativeSet(rRhs._derivativeSet),
      _buffer(rRhs._buffer.size()),
      _wrapped(rRhs._wrapped)
{
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(const AnsatzSpaceDerivative& rRhs) noexcept
requires hasStaticBasis
    : _ansatzSet(rRhs._ansatzSet),
      _derivativeSet(rRhs._derivativeSet),
      _buffer(),
      _wrapped(rRhs._wrapped)
{
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::evaluate(ConstSpan in, Span out) const {
    _wrapped.evaluate(in, out);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return _wrapped.size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::size() noexcept
requires hasStaticBasis {
    return View::size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (!hasStaticBasis) {
    return {_ansatzSet.data(), _ansatzSet.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const typename TScalarExpression::Derivative>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::derivativeSet() const noexcept
requires (!hasStaticBasis) {
    return {_derivativeSet.data(), _derivativeSet.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression,SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires hasStaticBasis {
    return {_ansatzSet.data(), _ansatzSet.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const typename TScalarExpression::Derivative,SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::derivativeSet() const noexcept
requires hasStaticBasis {
    return {_derivativeSet.data(), _derivativeSet.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::View
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (!hasStaticBasis) {
    return AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>(
        _ansatzSet,
        _derivativeSet);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::View
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires hasStaticBasis {
    return AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>(
        _ansatzSet,
        _derivativeSet);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView() noexcept
    : _set(),
      _buffer()
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(std::span<const TScalarExpression> ansatzSet) noexcept
requires (!hasStaticBasis)
    : _set(ansatzSet),
      _buffer()
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(std::span<const TScalarExpression,SetSize> ansatzSet) noexcept
requires hasStaticBasis
    : _set(ansatzSet),
      _buffer()
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(
    std::span<const TScalarExpression> ansatzSet,
    Span buffer)
requires (!hasStaticBasis)
    : _set(ansatzSet),
      _buffer() {
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(
    std::span<const TScalarExpression,SetSize> ansatzSet,
    std::span<Value,staticBufferSize> buffer) noexcept
requires hasStaticBasis
    : _set(ansatzSet),
      _buffer() {
    this->setBuffer(buffer);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceView<TScalarExpression,Dim,SetSize>::evaluate(ConstSpan in, Span out) const {
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


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return intPow(_set.size(), Dim);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::size() noexcept
requires hasStaticBasis {
    return intPow(SetSize, Dim);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::getMinBufferSize() const noexcept
requires (!hasStaticBasis) {
    // Size of the index buffer in bytes.
    constexpr unsigned indexBufferSize = Dim * sizeof(unsigned);

    // Offset of the value buffer from the begin of the buffer, in bytes.
    constexpr unsigned ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    // Required size of the value buffer in bytes.
    unsigned valueBufferSize = sizeof(Value) * intPow(_set.size(), Dim);

    return (ansatzBufferOffset + valueBufferSize) / sizeof(Value);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::getMinBufferSize() noexcept
requires hasStaticBasis {
    return staticBufferSize;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceView<TScalarExpression,Dim,SetSize>::setBuffer(Span buffer)
requires (!hasStaticBasis) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr void AnsatzSpaceView<TScalarExpression,Dim,SetSize>::setBuffer(std::span<Value,staticBufferSize> buffer) noexcept
requires hasStaticBasis {
    _buffer = buffer;
    auto indexBuffer = this->getIndexBuffer();
    std::fill_n(indexBuffer.data(), indexBuffer.size(), 0u);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<unsigned,Dim> AnsatzSpaceView<TScalarExpression,Dim,SetSize>::getIndexBuffer() const noexcept {
    return std::span<unsigned,Dim> {
        static_cast<unsigned*>(static_cast<void*>(_buffer.data())),
        Dim};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceView<TScalarExpression,Dim,SetSize>::Span
AnsatzSpaceView<TScalarExpression,Dim,SetSize>::getValueBuffer() const noexcept
requires (!hasStaticBasis) {
    return Span(
        _buffer.data() + ansatzBufferOffset / sizeof(Value),
        intPow(_set.size(), Dim));
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<
    typename AnsatzSpaceView<TScalarExpression,Dim,SetSize>::Value,
    AnsatzSpaceView<TScalarExpression,Dim,SetSize>::valueBufferSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::getValueBuffer() const noexcept
requires hasStaticBasis {
    return {
        _buffer.data() + ansatzBufferOffset / sizeof(Value),
        intPow(SetSize, Dim)};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression> AnsatzSpaceView<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (!hasStaticBasis) {
    return {_set.data(), _set.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<const TScalarExpression,SetSize> AnsatzSpaceView<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires hasStaticBasis {
    return {_set.data(), SetSize};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace() noexcept
    : AnsatzSpace(AnsatzSet {})
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSet&& rSet)
requires (!hasStaticBasis)
    : _set(std::move(rSet)),
      _buffer(),
      _wrapped()
{
    CIE_BEGIN_EXCEPTION_TRACING
    _wrapped = AnsatzSpaceView<TScalarExpression,Dim,SetSize>({_set.data(), _set.size()});
    _buffer.resize(_wrapped.getMinBufferSize());
    _wrapped.setBuffer({_buffer.data(), _buffer.size()});
    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSet&& rSet) noexcept
requires hasStaticBasis
    : _set(std::move(rSet)),
      _buffer(),
      _wrapped()
{
    _wrapped = AnsatzSpaceView<TScalarExpression,Dim,SetSize>(
        {_set.data(), _set.size()},
        {_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSet& rSet)
requires (!hasStaticBasis)
    : AnsatzSpace(AnsatzSet(rSet))
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSet& rSet) noexcept
requires hasStaticBasis
    : _set(std::move(rSet)),
      _buffer(),
      _wrapped()
{
    _wrapped = AnsatzSpaceView<TScalarExpression,Dim,SetSize>(
        {_set.data(), _set.size()},
        {_buffer.data(), _buffer.size()});
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSpace& rRhs)
requires (!hasStaticBasis)
    : AnsatzSpace(rRhs._set)
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSpace& rRhs) noexcept
requires hasStaticBasis
    : AnsatzSpace(rRhs._set)
{
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>&
AnsatzSpace<TScalarExpression,Dim,SetSize>::operator=(const AnsatzSpace& rRhs)
requires (!hasStaticBasis) {
    (*this) = AnsatzSpace(rRhs._set);
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>&
AnsatzSpace<TScalarExpression,Dim,SetSize>::operator=(const AnsatzSpace& rRhs) noexcept
requires hasStaticBasis {
    (*this) = AnsatzSpace(rRhs._set);
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpace<TScalarExpression,Dim,SetSize>::evaluate(ConstSpan in, Span out) const {
    _wrapped.evaluate(in, out);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpace<TScalarExpression,Dim,SetSize>::Derivative
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeDerivative() const
requires (!hasStaticBasis) {
    return Derivative(this->ansatzSet());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr typename AnsatzSpace<TScalarExpression,Dim,SetSize>::Derivative
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeDerivative() const noexcept
requires hasStaticBasis {
    return Derivative(this->ansatzSet());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return _wrapped.size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::size() noexcept
requires hasStaticBasis {
    return View::size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression> AnsatzSpace<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (!hasStaticBasis) {
    return {_set.data(), _set.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<const TScalarExpression,SetSize> AnsatzSpace<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires hasStaticBasis {
    return {_set.data(), _set.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpace<TScalarExpression,Dim,SetSize>::View
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (!hasStaticBasis) {
    return AnsatzSpaceView<TScalarExpression,Dim,SetSize>(_set);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr typename AnsatzSpace<TScalarExpression,Dim,SetSize>::View
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires hasStaticBasis {
    return AnsatzSpaceView<TScalarExpression,Dim,SetSize>(_set);
}


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::header(Ref<XMLElement> rElement)
{
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement defaultData = rElement.addChild("default");
    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::operator()(Ref<XMLElement> rElement,
                                                                                Ref<const maths::AnsatzSpace<TScalarExpression,Dim,SetSize>> rInstance)
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


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onElementBegin(Ptr<void> pThis,
                                                                                      std::string_view elementName,
                                                                                      [[maybe_unused]] std::span<GraphML::AttributePair> attributes) {
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


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onText(Ptr<void>,
                                                                              std::string_view) {
    CIE_THROW(
        Exception,
        "Unexpected text block while parsing AnsatzSpace in GraphML."
    )
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onElementEnd(Ptr<void> pThis,
                                                                                            std::string_view elementName) {
    CIE_BEGIN_EXCEPTION_TRACING
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = Value(std::move(rThis._set));
    rThis.template release<Deserializer>(&rThis, elementName);
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem::io
