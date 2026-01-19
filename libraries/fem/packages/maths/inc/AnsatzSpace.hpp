#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/maths/inc/power.hpp"

// --- STL Includes ---
#include <array>


namespace cie::fem::maths {


template <class, unsigned, std::size_t>
class AnsatzSpace;


template <class, unsigned, std::size_t>
class AnsatzSpaceDerivative;


template <class TScalarExpression, unsigned Dim, std::size_t SetSize = 0ul>
class AnsatzSpaceDerivativeView : public ExpressionTraits<typename TScalarExpression::Value> {
public:
    constexpr static inline bool hasStaticBasis = (SetSize != 0ul);

    constexpr static inline std::size_t BasisCount = SetSize;

    static constexpr unsigned Dimension = Dim;

    using typename ExpressionTraits<typename TScalarExpression::Value>::Value;

    using typename ExpressionTraits<Value>::ConstSpan;

    using typename ExpressionTraits<Value>::Span;

private:
    static inline constexpr std::size_t indexBufferSize = Dimension * sizeof(unsigned);

    static inline constexpr std::size_t ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    static inline constexpr std::size_t valueBufferSize = intPow(BasisCount, Dim) * sizeof(Value);

    static inline constexpr std::size_t derivativeBufferOffset = ansatzBufferOffset + valueBufferSize;

public:
    static inline constexpr std::size_t staticBufferSize = (derivativeBufferOffset + valueBufferSize) / sizeof(Value);

    constexpr AnsatzSpaceDerivativeView() noexcept = default;

    constexpr AnsatzSpaceDerivativeView(std::span<const TScalarExpression,SetSize> ansatzSet,
                                        std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet) noexcept
    requires hasStaticBasis;

    constexpr AnsatzSpaceDerivativeView(std::span<const TScalarExpression,SetSize> ansatzSet,
                                        std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet,
                                        std::span<Value,staticBufferSize> buffer) noexcept
    requires hasStaticBasis;

    AnsatzSpaceDerivativeView(std::span<const TScalarExpression> ansatzSet,
                              std::span<const typename TScalarExpression::Derivative> derivativeSet)
    requires (!hasStaticBasis);

    AnsatzSpaceDerivativeView(std::span<const TScalarExpression> ansatzSet,
                              std::span<const typename TScalarExpression::Derivative> derivativeSet,
                              Span buffer)
    requires (!hasStaticBasis);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept
    requires (!hasStaticBasis);

    unsigned getMinBufferSize() const noexcept
    requires (!hasStaticBasis);

    void setBuffer(Span buffer)
    requires (!hasStaticBasis);

    static constexpr unsigned size() noexcept
    requires hasStaticBasis;

    static constexpr unsigned getMinBufferSize() noexcept
    requires hasStaticBasis;

    constexpr void setBuffer(std::span<Value,staticBufferSize> buffer) noexcept
    requires hasStaticBasis;

    constexpr std::span<const TScalarExpression> ansatzSet() const noexcept;

    constexpr std::span<const typename TScalarExpression::Derivative> derivativeSet() const noexcept;

private:
    friend class AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>;

    constexpr std::span<unsigned,Dim> getIndexBuffer() const noexcept;

    Span getAnsatzBuffer() const noexcept
    requires (!hasStaticBasis);

    Span getDerivativeBuffer() const noexcept
    requires (!hasStaticBasis);

    constexpr std::span<Value,valueBufferSize> getAnsatzBuffer() const noexcept
    requires hasStaticBasis;

    constexpr std::span<Value,valueBufferSize> getDerivativeBuffer() const noexcept
    requires hasStaticBasis;

    std::conditional_t<
        hasStaticBasis,
        std::span<const TScalarExpression,SetSize>,
        std::span<const TScalarExpression>
    > _ansatzSet;

    std::conditional_t<
        hasStaticBasis,
        std::span<const typename TScalarExpression::Derivative,SetSize>,
        std::span<const typename TScalarExpression::Derivative>
    > _derivativeSet;

    std::conditional_t<
        hasStaticBasis,
        std::span<Value,staticBufferSize>,
        Span
    > _buffer;
}; // class AnsatzSpaceDerivativeView


template <class TScalarExpression, unsigned Dim, std::size_t SetSize = 0ul>
class AnsatzSpaceDerivative : public ExpressionTraits<typename TScalarExpression::Value> {
public:
    using View = AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize>;

    static inline constexpr bool hasStaticBasis = View::hasStaticBasis;

    static inline constexpr std::size_t BasisCount = SetSize;

    static constexpr unsigned Dimension = Dim;

    using typename ExpressionTraits<typename TScalarExpression::Value>::Value;

    using typename ExpressionTraits<Value>::ConstSpan;

    using typename ExpressionTraits<Value>::Span;

    constexpr AnsatzSpaceDerivative() noexcept = default;

    AnsatzSpaceDerivative(std::span<const TScalarExpression> ansatzSet)
    requires (!hasStaticBasis);

    constexpr AnsatzSpaceDerivative(std::span<const TScalarExpression,SetSize> ansatzSet) noexcept
    requires hasStaticBasis;

    AnsatzSpaceDerivative(AnsatzSpaceDerivative&&) noexcept
    requires (!hasStaticBasis) = default;

    constexpr AnsatzSpaceDerivative(AnsatzSpaceDerivative&&) noexcept
    requires hasStaticBasis = default;

    AnsatzSpaceDerivative(const AnsatzSpaceDerivative&)
    requires (!hasStaticBasis);

    constexpr AnsatzSpaceDerivative(const AnsatzSpaceDerivative&) noexcept
    requires hasStaticBasis;

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept
    requires (!hasStaticBasis);

    View makeView() const noexcept
    requires (!hasStaticBasis);

    std::span<const TScalarExpression> ansatzSet() const noexcept
    requires (!hasStaticBasis);

    std::span<const typename TScalarExpression::Derivative> derivativeSet() const noexcept
    requires (!hasStaticBasis);

    constexpr static unsigned size() noexcept
    requires hasStaticBasis;

    constexpr View makeView() const noexcept
    requires hasStaticBasis;

    constexpr std::span<const TScalarExpression,SetSize> ansatzSet() const noexcept
    requires hasStaticBasis;

    constexpr std::span<const typename TScalarExpression::Derivative,SetSize> derivativeSet() const noexcept
    requires hasStaticBasis;

private:
    friend class AnsatzSpace<TScalarExpression,Dim,SetSize>;

    std::conditional_t<
        hasStaticBasis,
        std::array<TScalarExpression,BasisCount>,
        DynamicArray<TScalarExpression>
    > _ansatzSet;

    std::conditional_t<
        hasStaticBasis,
        std::array<typename TScalarExpression::Derivative,BasisCount>,
        DynamicArray<typename TScalarExpression::Derivative>
    > _derivativeSet;

    std::conditional_t<
        hasStaticBasis,
        std::array<Value,View::staticBufferSize>,
        DynamicArray<Value>
    > _buffer;

    AnsatzSpaceDerivativeView<TScalarExpression,Dim,SetSize> _wrapped;
}; // class AnsatzSpaceDerivative


/** @brief A set of multidimensional functions constructed from the cartesian product of a set of scalar basis functions.
 */
template <class TScalarExpression, unsigned Dim, std::size_t SetSize = 0ul>
class AnsatzSpaceView : public ExpressionTraits<typename TScalarExpression::Value> {
public:
    using Base = ExpressionTraits<typename TScalarExpression::Value>;

    constexpr static inline bool hasStaticBasis = (SetSize != 0ul);

    constexpr static inline std::size_t BasisCount = SetSize;

    static constexpr unsigned Dimension = Dim;

    using typename Base::Value;

    using typename Base::Span;

    using typename Base::ConstSpan;

private:
    static inline constexpr std::size_t indexBufferSize = Dim * sizeof(unsigned);

    static inline constexpr std::size_t ansatzBufferOffset = (indexBufferSize / sizeof(Value) + (indexBufferSize % sizeof(Value) != 0)) * sizeof(Value);

    static inline constexpr std::size_t valueBufferSize = intPow(BasisCount, Dim) * sizeof(Value);

public:
    static inline constexpr std::size_t staticBufferSize = valueBufferSize / sizeof(Value);

    constexpr AnsatzSpaceView() noexcept;

    constexpr AnsatzSpaceView(std::span<const TScalarExpression> ansatzSet) noexcept
    requires (!hasStaticBasis);

    constexpr AnsatzSpaceView(std::span<const TScalarExpression,SetSize> ansatzSet) noexcept
    requires hasStaticBasis;

    AnsatzSpaceView(std::span<const TScalarExpression> ansatzSet,
                    Span buffer)
    requires (!hasStaticBasis);

    constexpr AnsatzSpaceView(std::span<const TScalarExpression,SetSize> ansatzSet,
                              std::span<Value,staticBufferSize> buffer) noexcept
    requires hasStaticBasis;

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept
    requires (!hasStaticBasis);

    unsigned getMinBufferSize() const noexcept
    requires (!hasStaticBasis);

    void setBuffer(Span buffer)
    requires (!hasStaticBasis);

    std::span<const TScalarExpression> ansatzSet() const noexcept
    requires (!hasStaticBasis);

    static constexpr unsigned size() noexcept
    requires hasStaticBasis;

    static constexpr unsigned getMinBufferSize() noexcept
    requires hasStaticBasis;

    constexpr void setBuffer(std::span<Value,staticBufferSize> buffer) noexcept
    requires hasStaticBasis;

    constexpr std::span<const TScalarExpression,SetSize> ansatzSet() const noexcept
    requires hasStaticBasis;

private:
    constexpr std::span<unsigned,Dim> getIndexBuffer() const noexcept;

    Span getValueBuffer() const noexcept
    requires (!hasStaticBasis);

    constexpr std::span<Value,valueBufferSize> getValueBuffer() const noexcept
    requires hasStaticBasis;

    std::conditional_t<
        hasStaticBasis,
        std::span<const TScalarExpression,SetSize>,
        std::span<const TScalarExpression>
    > _set;

    std::conditional_t<
        hasStaticBasis,
        std::span<Value,staticBufferSize>,
        Span
    > _buffer;
}; // class AnsatzSpaceView



/** @brief A set of multidimensional functions constructed from the cartesian product of a set of scalar basis functions.
 */
template <class TScalarExpression, unsigned Dim, std::size_t SetSize = 0ul>
class AnsatzSpace : public ExpressionTraits<typename TScalarExpression::Value> {
public:
    using Base = ExpressionTraits<typename TScalarExpression::Value>;

    using View = AnsatzSpaceView<TScalarExpression,Dim,SetSize>;

    static inline constexpr std::size_t BasisCount = View::BasisCount;

    static inline constexpr bool hasStaticBasis = View::hasStaticBasis;

    static constexpr unsigned Dimension = Dim;

    using typename Base::Value;

    using typename Base::Span;

    using typename Base::ConstSpan;

    using AnsatzSet = std::conditional_t<
        hasStaticBasis,
        std::array<TScalarExpression,BasisCount>,
        DynamicArray<TScalarExpression>>;

    using Derivative = AnsatzSpaceDerivative<TScalarExpression,Dim,SetSize>;

public:
    constexpr AnsatzSpace() noexcept;

    AnsatzSpace(AnsatzSet&& rSet)
    requires (!hasStaticBasis);

    constexpr AnsatzSpace(AnsatzSet&& rSet) noexcept
    requires hasStaticBasis;

    AnsatzSpace(const AnsatzSet& rSet)
    requires (!hasStaticBasis);

    constexpr AnsatzSpace(const AnsatzSet& rSet) noexcept
    requires hasStaticBasis;

    AnsatzSpace(AnsatzSpace&&) noexcept
    requires (!hasStaticBasis) = default;

    constexpr AnsatzSpace(AnsatzSpace&&) noexcept
    requires hasStaticBasis = default;

    AnsatzSpace(const AnsatzSpace& rRhs)
    requires (!hasStaticBasis);

    constexpr AnsatzSpace(const AnsatzSpace& rRhs) noexcept
    requires hasStaticBasis;

    AnsatzSpace& operator=(AnsatzSpace&&) noexcept = default;

    AnsatzSpace& operator=(const AnsatzSpace& rRhs)
    requires (!hasStaticBasis);

    constexpr AnsatzSpace& operator=(const AnsatzSpace& rRhs) noexcept
    requires hasStaticBasis;

    void evaluate(ConstSpan in, Span out) const;

    Derivative makeDerivative() const
    requires (!hasStaticBasis);

    constexpr Derivative makeDerivative() const noexcept
    requires hasStaticBasis;

    unsigned size() const noexcept
    requires (!hasStaticBasis);

    constexpr static unsigned size() noexcept
    requires hasStaticBasis;

    std::span<const TScalarExpression> ansatzSet() const noexcept
    requires (!hasStaticBasis);

    constexpr std::span<const TScalarExpression,SetSize> ansatzSet() const noexcept
    requires hasStaticBasis;

    View makeView() const noexcept
    requires (!hasStaticBasis);

    constexpr View makeView() const noexcept
    requires hasStaticBasis;

private:
    friend class AnsatzSpaceDerivative<TScalarExpression,Dim>;

    AnsatzSet _set;

    std::conditional_t<
        hasStaticBasis,
        std::array<Value,View::staticBufferSize>,
        DynamicArray<Value>
    > _buffer;

    AnsatzSpaceView<TScalarExpression,Dim,SetSize> _wrapped;
}; // class AnsatzSpace


} // namespace cie::fem::maths



// --- IO --- //



namespace cie::fem::io {


template <class TScalarExpression, unsigned Dimension, std::size_t SetSize>
struct io::GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>>
{
    using Value = maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>;

    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>> rInstance);
}; // struct GraphML::Serializer<AnsatzSpace>


template <class TScalarExpression, unsigned Dimension, std::size_t SetSize>
struct io::GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>>
    : public io::GraphML::DeserializerBase<maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>>
{
    using Value = maths::AnsatzSpace<TScalarExpression,Dimension,SetSize>;

    using io::GraphML::DeserializerBase<Value>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);

private:
    typename Value::AnsatzSet _set;
}; // struct GraphML::Deserializer<AnsatzSpace>


} // namespace cie::fem::io

#include "packages/maths/impl/AnsatzSpace_impl.hpp"
