#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>
#include <numeric>


namespace cie::fem::maths {


namespace impl {

/// @brief A struct defining the buffer layout of the dynamic version of @ref AnsatzSpace.
/// @details Assuming the ansatz set consists of @p n functions, the buffer layout of @ref AnsatzSpace looks like this:
///          @code
///          | _ _ _ | _ _ _ |
///              ^       ^
///              |       |
///              |       + nested buffer required by the ansatz set (max)
///              + buffer for the values of the basis functions (@p n scalars)
///          @endcode
template <class TScalarExpression>
struct AnsatzTraits {
    using Value = typename TScalarExpression::Value;

    static constexpr std::size_t ansatzBufferOffset([[maybe_unused]] std::span<const TScalarExpression> set) noexcept {
        return 0u;
    }

    static constexpr std::size_t ansatzBufferSize(std::span<const TScalarExpression> set) noexcept {
        return std::accumulate(
            set.begin(),
            set.end(),
            0u,
            [] (unsigned left, Ref<const TScalarExpression> f) {return left + f.size();});
    }

    static constexpr std::size_t nestedBufferOffset(std::span<const TScalarExpression> set) noexcept {
        return AnsatzTraits::ansatzBufferOffset(set) + AnsatzTraits::ansatzBufferSize(set);
    }

    static constexpr std::size_t nestedBufferSize(std::span<const TScalarExpression> set) noexcept {
        return set.empty()
            ? 0ul
            : std::max_element(
                set.begin(),
                set.end(),
                [] (Ref<const TScalarExpression> f) {return f.bufferSize();}
              )->bufferSize();
    }

    static constexpr std::size_t bufferSize(std::span<const TScalarExpression> set) noexcept {
        return AnsatzTraits::nestedBufferOffset(set) + AnsatzTraits::nestedBufferSize(set);
    }

    static constexpr typename ExpressionTraits<Value>::Span getAnsatzBuffer(
        std::span<const TScalarExpression> set,
        typename ExpressionTraits<Value>::BufferView buffer) noexcept {
            return typename ExpressionTraits<Value>::Span(
                buffer.data() + AnsatzTraits::ansatzBufferOffset(set),
                AnsatzTraits::ansatzBufferSize(set));
    }

    static constexpr typename ExpressionTraits<Value>::Span getNestedBuffer(
        std::span<const TScalarExpression> set,
        typename ExpressionTraits<Value>::BufferView buffer) noexcept {
        return typename ExpressionTraits<Value>::Span(
            buffer.data() + AnsatzTraits::nestedBufferOffset(set),
            AnsatzTraits::nestedBufferSize(set));
    }
}; // struct AnsatzTraits


/// @brief A struct defining the buffer layout of the dynamic version of @ref AnsatzSpaceDerivative.
/// @details Assuming the ansatz set consists of @p n functions, the buffer layout of @ref AnsatzSpaceDerivative looks like this:
///          @code
///          | _ _ _ | _ _ _ | _ _ _ |
///              ^       ^       ^
///              |       |       |
///              |       |       + nested buffer required by the ansatz set and its derivative (max)
///              |       + buffer for the values of ansatz
///              + buffer for the values of the basis functions (@p n scalars)
///          @endcode
template <class TScalarExpression>
struct AnsatzDerivativeTraits {
    using Value = typename TScalarExpression::Value;

    static constexpr std::size_t ansatzBufferOffset(
        [[maybe_unused]] std::span<const TScalarExpression> set,
        [[maybe_unused]] std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
            return 0u;
    }

    static constexpr std::size_t ansatzBufferSize(
        std::span<const TScalarExpression> set,
        [[maybe_unused]] std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
        return std::accumulate(
            set.begin(),
            set.end(),
            0u,
            [] (unsigned left, Ref<const TScalarExpression> f) {return left + f.size();});
    }

    static constexpr std::size_t derivativeBufferOffset(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
            return AnsatzDerivativeTraits::ansatzBufferOffset(set, derivativeSet)
                 + AnsatzDerivativeTraits::ansatzBufferSize(set, derivativeSet);
    }

    static constexpr std::size_t derivativeBufferSize(
        [[maybe_unused]] std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
        return std::accumulate(
            derivativeSet.begin(),
            derivativeSet.end(),
            0u,
            [] (unsigned left, Ref<const typename TScalarExpression::Derivative> f) {return left + f.size();});
    }

    static constexpr std::size_t nestedBufferOffset(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
            return AnsatzDerivativeTraits::derivativeBufferOffset(set, derivativeSet)
                 + AnsatzDerivativeTraits::derivativeBufferSize(set, derivativeSet);
    }

    static constexpr std::size_t nestedBufferSize(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
        return std::max<unsigned>(
            set.empty()
                ? 0u
                : std::max_element(
                    set.begin(),
                    set.end(),
                    [] (Ref<const TScalarExpression> f) {return f.bufferSize();}
                )->bufferSize(),
            derivativeSet.empty()
                ? 0u
                : std::max_element(
                    derivativeSet.begin(),
                    derivativeSet.end(),
                    [] (Ref<const TScalarExpression> f) {return f.bufferSize();}
                )->bufferSize()
            );
    }

    static constexpr std::size_t bufferSize(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet) noexcept {
            return AnsatzDerivativeTraits::nestedBufferOffset(set, derivativeSet)
                 + AnsatzDerivativeTraits::nestedBufferSize(set, derivativeSet);
    }

    static constexpr typename ExpressionTraits<Value>::Span getAnsatzBuffer(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet,
        typename ExpressionTraits<Value>::BufferView buffer) noexcept {
            return typename ExpressionTraits<Value>::Span(
                buffer.data() + AnsatzDerivativeTraits::ansatzBufferOffset(set, derivativeSet),
                AnsatzDerivativeTraits::ansatzBufferSize(set, derivativeSet));
    }

    static constexpr typename ExpressionTraits<Value>::Span getDerivativeBuffer(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet,
        typename ExpressionTraits<Value>::BufferView buffer) noexcept {
            return typename ExpressionTraits<Value>::Span(
                buffer.data() + AnsatzDerivativeTraits::derivativeBufferOffset(set, derivativeSet),
                AnsatzDerivativeTraits::derivativeBufferSize(set, derivativeSet));
    }

    static constexpr typename ExpressionTraits<Value>::Span getNestedBuffer(
        std::span<const TScalarExpression> set,
        std::span<const typename TScalarExpression::Derivative> derivativeSet,
        typename ExpressionTraits<Value>::BufferView buffer) noexcept {
            return typename ExpressionTraits<Value>::Span(
                buffer.data() + AnsatzDerivativeTraits::nestedBufferOffset(set, derivativeSet),
                AnsatzDerivativeTraits::nestedBufferSize(set, derivativeSet));
    }
}; // struct AnsatzTraits

} // namespace impl


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView() noexcept
    : _ansatzSet(static_cast<const TScalarExpression*>(nullptr), SetSize),
      _derivativeSet(static_cast<const typename TScalarExpression::Derivative*>(nullptr), SetSize)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression,SetSize> ansatzSet,
    std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet) noexcept
requires (hasStaticBasis)
    : _ansatzSet(ansatzSet),
      _derivativeSet(derivativeSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivativeView(
    std::span<const TScalarExpression> ansatzSet,
    std::span<const typename TScalarExpression::Derivative> derivativeSet)
requires (!hasStaticBasis)
    : _ansatzSet(ansatzSet),
      _derivativeSet(derivativeSet) {
        assert(_derivativeSet.size() == _ansatzSet.size());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        const unsigned setSize = _ansatzSet.size();

        // Santity checks.
        assert(in.size() == Dim); /// @todo pass an ansatz mask here
        assert(setSize == _derivativeSet.size());
        assert(out.size() == this->size());

        // Construct a buffer for basis indices.
        std::array<std::uint16_t,Dim> indexBuffer;
        std::fill(indexBuffer.begin(), indexBuffer.end(), static_cast<std::uint16_t>(0));

        // Declare required buffers that will be filled set
        // depending on whether the ansatz set is dynamically
        // sized.
        Span ansatzBuffer, derivativeBuffer;
        BufferSpan nestedBuffer;

        constexpr std::size_t staticNestedBufferSize = std::max<std::size_t>(
            StaticExpressionSize<TScalarExpression>::bufferSize,
            StaticExpressionSize<typename TScalarExpression::Derivative>::bufferSize);
        [[maybe_unused]] std::array<typename TScalarExpression::Value,SetSize> ansatzBufferArray, derivativeBufferArray;
        [[maybe_unused]] std::array<typename BufferSpan::value_type,staticNestedBufferSize> nestedBufferArray;

        // Set the required buffers.
        if constexpr (hasStaticBasis) {
            ansatzBuffer        = ansatzBufferArray;
            derivativeBuffer    = derivativeBufferArray;
            nestedBuffer        = nestedBufferArray;
        } else {
            ansatzBuffer = impl::AnsatzDerivativeTraits<TScalarExpression>::getAnsatzBuffer(
                _ansatzSet,
                _derivativeSet,
                buffer);
            derivativeBuffer = impl::AnsatzDerivativeTraits<TScalarExpression>::getDerivativeBuffer(
                _ansatzSet,
                _derivativeSet,
                buffer);
            nestedBuffer = impl::AnsatzDerivativeTraits<TScalarExpression>::getNestedBuffer(
                _ansatzSet,
                _derivativeSet,
                buffer);
        }

        // Fill the value and derivative buffers
        {
            Ptr<Value> pValue      = ansatzBuffer.data();
            Ptr<Value> pDerivative = derivativeBuffer.data();

            for (auto c : in) {
                ConstSpan scalarIn {&c, 1};
                for (const auto& rScalarExpression : _ansatzSet) {
                    rScalarExpression.evaluate(
                        scalarIn,
                        {pValue, 1},
                        {nestedBuffer.data(), rScalarExpression.size()});
                    ++pValue;
                } // for rScalarExpression in ansatzSet
                for (const auto& rScalarExpression : _derivativeSet) {
                    rScalarExpression.evaluate(
                        scalarIn,
                        {pDerivative, 1},
                        {nestedBuffer.data(), rScalarExpression.size()});
                    ++pDerivative;
                } // for rScalarExpression in derivativeSet
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
requires (hasStaticBasis) {
    return intPow(SetSize, Dim) * Dim;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::bufferSize() const noexcept
requires (!hasStaticBasis) {
    return impl::AnsatzDerivativeTraits<TScalarExpression>::bufferSize(_ansatzSet, _derivativeSet);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>::bufferSize() noexcept
requires (hasStaticBasis) {
    return std::max<unsigned>(
        TScalarExpression::bufferSize(),
        TScalarExpression::Derivative::bufferSize());
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
requires (hasStaticBasis) {
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
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(
    AnsatzSpaceDerivative&& rRhs) noexcept
requires (!hasStaticBasis)
    : _ansatzSet(std::move(rRhs._ansatzSet)),
      _derivativeSet(std::move(rRhs._derivativeSet))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(
    AnsatzSpaceDerivative&& rRhs) noexcept
requires (hasStaticBasis)
    : _ansatzSet(rRhs._ansatzSet),
      _derivativeSet(rRhs._derivativeSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet)
requires (!hasStaticBasis)
    : _ansatzSet(ansatzSet.begin(), ansatzSet.end()),
      _derivativeSet(ansatzSet.size())
{
    std::transform(
        _ansatzSet.begin(),
        _ansatzSet.end(),
        _derivativeSet.begin(),
        [](Ref<const TScalarExpression> rAnsatzFunction){
             return rAnsatzFunction.makeDerivative();
        });
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(const AnsatzSpaceDerivative& rRhs)
requires (!hasStaticBasis)
    : _ansatzSet(rRhs._ansatzSet),
      _derivativeSet(rRhs._derivativeSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::AnsatzSpaceDerivative(const AnsatzSpaceDerivative& rRhs) noexcept
requires (hasStaticBasis)
    : _ansatzSet(rRhs._ansatzSet),
      _derivativeSet(rRhs._derivativeSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>&
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::operator=(AnsatzSpaceDerivative&& rRhs) noexcept
requires (!hasStaticBasis) {
    _ansatzSet = std::move(rRhs._ansatzSet);
    _derivativeSet = std::move(rRhs._derivativeSet);
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>&
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::operator=(AnsatzSpaceDerivative&& rRhs) noexcept
requires (hasStaticBasis) {
    _ansatzSet = rRhs._ansatzSet;
    _derivativeSet = rRhs._derivativeSet;
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>&
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::operator=(const AnsatzSpaceDerivative& rRhs)
requires (!hasStaticBasis) {
    _ansatzSet = rRhs._ansatzSet;
    _derivativeSet = rRhs._derivativeSet;
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>&
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::operator=(const AnsatzSpaceDerivative& rRhs) noexcept
requires (hasStaticBasis) {
    _ansatzSet = rRhs._ansatzSet;
    _derivativeSet = rRhs._derivativeSet;
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        return this->makeView().evaluate(
            in,
            out,
            buffer);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return this->makeView().size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::size() noexcept
requires (hasStaticBasis) {
    return View::size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::bufferSize() const noexcept
requires (!hasStaticBasis) {
    return this->makeView().bufferSize();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::bufferSize() noexcept
requires (hasStaticBasis) {
    return View::bufferSize();
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
requires (hasStaticBasis) {
    return std::span<const TScalarExpression,SetSize>(
        _ansatzSet.data(),
        _ansatzSet.size());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const typename TScalarExpression::Derivative,SetSize>
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::derivativeSet() const noexcept
requires (hasStaticBasis) {
    return std::span<const typename TScalarExpression::Derivative,SetSize>(
        _derivativeSet.data(),
        _derivativeSet.size());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::View
AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (!hasStaticBasis) {
    return View(
        _ansatzSet,
        _derivativeSet);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::View
constexpr AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (hasStaticBasis) {
    return View(
        _ansatzSet,
        _derivativeSet);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView() noexcept
    : _set(static_cast<const TScalarExpression*>(nullptr), SetSize)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(std::span<const TScalarExpression> ansatzSet) noexcept
requires (!hasStaticBasis)
    : _set(ansatzSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpaceView<TScalarExpression,Dim,SetSize>::AnsatzSpaceView(std::span<const TScalarExpression,SetSize> ansatzSet) noexcept
requires (hasStaticBasis)
    : _set(ansatzSet)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void AnsatzSpaceView<TScalarExpression,Dim,SetSize>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        // Sanity checks
        assert(in.size() == Dim);
        assert(out.size() == this->size());

        // Construct a buffer for basis indices.
        std::array<std::uint16_t,Dim> indexBuffer;
        std::fill(indexBuffer.begin(), indexBuffer.end(), static_cast<std::uint16_t>(0));

        // Declare required buffers that will be filled set
        // depending on whether the ansatz set is dynamically
        // sized.
        Span ansatzBuffer;
        BufferSpan nestedBuffer;

        constexpr std::size_t staticNestedBufferSize = StaticExpressionSize<TScalarExpression>::bufferSize;
        [[maybe_unused]] std::array<typename TScalarExpression::Value,SetSize> ansatzBufferArray;
        [[maybe_unused]] std::array<typename BufferSpan::value_type,staticNestedBufferSize> nestedBufferArray;

        // Set the required buffers.
        if constexpr (hasStaticBasis) {
            ansatzBuffer = ansatzBufferArray;
            nestedBuffer = nestedBufferArray;
        } else {
            ansatzBuffer = impl::AnsatzTraits<TScalarExpression>::getAnsatzBuffer(
                _set,
                buffer);
            nestedBuffer = impl::AnsatzTraits<TScalarExpression>::getNestedBuffer(
                _set,
                buffer);
        }

        // Fill the value buffer
        auto pValue = ansatzBuffer.data();
        for (unsigned iDim=0u; iDim<Dim; ++iDim) {
            for (const auto& rScalarExpression : _set) {
                rScalarExpression.evaluate(
                    {in.data() + iDim, 1},
                    {pValue++, 1},
                    nestedBuffer);
            } // for rScalarExpression in _set
        } // for iDim in range(Dim)

        const unsigned setSize = _set.size();
        auto itOut = out.data();
        do {
            // Compute product of bases
            *itOut = static_cast<Value>(1);
            for (unsigned iIndex=0; iIndex<indexBuffer.size(); ++iIndex) {
                *itOut *= ansatzBuffer[indexBuffer[iIndex] + iIndex * setSize];
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
requires (hasStaticBasis) {
    return intPow(SetSize, Dim);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::bufferSize() const noexcept
requires (!hasStaticBasis) {
    return impl::AnsatzTraits<TScalarExpression>::bufferSize(_set);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpaceView<TScalarExpression,Dim,SetSize>::bufferSize() noexcept
requires (hasStaticBasis) {
    return TScalarExpression::bufferSize();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression> AnsatzSpaceView<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (!hasStaticBasis) {
    return {_set.data(), _set.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<const TScalarExpression,SetSize> AnsatzSpaceView<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (hasStaticBasis) {
    return {_set.data(), SetSize};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace() noexcept
    : AnsatzSpace(AnsatzSet {})
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSet&& rSet)
requires (!hasStaticBasis)
    : _set(std::move(rSet))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSet&& rSet) noexcept
requires (hasStaticBasis)
    : _set(std::move(rSet))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSet& rSet)
requires (!hasStaticBasis)
    : AnsatzSpace(AnsatzSet(rSet))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSet& rSet) noexcept
requires (hasStaticBasis)
    : _set(std::move(rSet))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSpace&& rRhs) noexcept
requires (!hasStaticBasis)
    : AnsatzSpace(std::move(rRhs._set))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(AnsatzSpace&& rRhs) noexcept
requires (hasStaticBasis)
    : AnsatzSpace(std::move(rRhs._set))
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSpace& rRhs)
requires (!hasStaticBasis)
    : AnsatzSpace(rRhs._set)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>::AnsatzSpace(const AnsatzSpace& rRhs) noexcept
requires (hasStaticBasis)
    : AnsatzSpace(rRhs._set)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
AnsatzSpace<TScalarExpression,Dim,SetSize>&
AnsatzSpace<TScalarExpression,Dim,SetSize>::operator=(AnsatzSpace&& rRhs) noexcept
requires (!hasStaticBasis) {
    _set = std::move(rRhs._set);
    return *this;
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr AnsatzSpace<TScalarExpression,Dim,SetSize>&
AnsatzSpace<TScalarExpression,Dim,SetSize>::operator=(AnsatzSpace&& rRhs) noexcept
requires (hasStaticBasis) {
    _set = rRhs._set;
    return *this;
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
requires (hasStaticBasis) {
    (*this) = AnsatzSpace(rRhs._set);
    return *this;
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
requires (hasStaticBasis) {
    return Derivative(this->ansatzSet());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::size() const noexcept
requires (!hasStaticBasis) {
    return this->makeView().size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::size() noexcept
requires (hasStaticBasis) {
    return View::size();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::bufferSize() const noexcept
requires (!hasStaticBasis) {
    return this->makeView().bufferSize();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr unsigned AnsatzSpace<TScalarExpression,Dim,SetSize>::bufferSize() noexcept
requires (hasStaticBasis) {
    return View::bufferSize();
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
std::span<const TScalarExpression> AnsatzSpace<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (!hasStaticBasis) {
    return {_set.data(), _set.size()};
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr std::span<const TScalarExpression,SetSize> AnsatzSpace<TScalarExpression,Dim,SetSize>::ansatzSet() const noexcept
requires (hasStaticBasis) {
    return std::span<const TScalarExpression,SetSize>(_set.data(), _set.size());
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
typename AnsatzSpace<TScalarExpression,Dim,SetSize>::View
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (!hasStaticBasis) {
    return View(_set);
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
constexpr typename AnsatzSpace<TScalarExpression,Dim,SetSize>::View
AnsatzSpace<TScalarExpression,Dim,SetSize>::makeView() const noexcept
requires (hasStaticBasis) {
    return View(_set);
}


} // namespace cie::fem::maths


namespace cie::fem::io {


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::header(Ref<XMLElement> rElement) {
    CIE_BEGIN_EXCEPTION_TRACING
        GraphML::XMLElement defaultData = rElement.addChild("default");
    CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::operator()(
    Ref<XMLElement> rElement,
    Ref<const maths::AnsatzSpace<TScalarExpression,Dim,SetSize>> rInstance) {
        CIE_BEGIN_EXCEPTION_TRACING
            using SubSerializer = GraphML::Serializer<std::span<const TScalarExpression>>;
            SubSerializer subSerializer;
            GraphML::XMLElement subElement = rElement.addChild("ansatz-space");
            subSerializer(subElement, rInstance.ansatzSet());
        CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onElementBegin(
    Ptr<void> pThis,
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
                SubDeserializer::onElementEnd});
        CIE_END_EXCEPTION_TRACING
}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onText(
    Ptr<void>,
    std::string_view)
{}


template <class TScalarExpression, unsigned Dim, std::size_t SetSize>
void GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dim,SetSize>>::onElementEnd(
    Ptr<void> pThis,
    std::string_view elementName) {
        CIE_BEGIN_EXCEPTION_TRACING
            Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
            rThis.instance() = Value(std::move(rThis._set));
            rThis.template release<Deserializer>(&rThis, elementName);
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem::io
