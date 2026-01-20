#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
class TransformedIntegrand : public maths::ExpressionTraits<typename TIntegrand::Value> {
public:
    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::ConstSpan;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Span;

    using Integrand = TIntegrand;

    using Jacobian = TJacobian;

public:
    constexpr static inline bool isBuffered =
        maths::BufferedExpression<TIntegrand>
        || maths::BufferedExpression<TJacobian>;

    TransformedIntegrand() noexcept {}

    TransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         RightRef<TJacobian> rInverseJacobian) noexcept
        : _integrand(std::move(rIntegrand)),
          _inverseJacobian(std::move(rInverseJacobian))
    {}

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TIntegrand>) {
        return _integrand.size();}

    static constexpr unsigned size() noexcept
    requires (maths::StaticExpression<TIntegrand>) {
        return TIntegrand::size();}

    void evaluate(ConstSpan in, Span out) const {
        const Value determinant = std::abs(_inverseJacobian.evaluateDeterminant(in));
        const Value scale = static_cast<Value>(1) / determinant;
        _integrand.evaluate(in, out);
        for (Value& rComponent : out) {
            rComponent *= scale;
        }
    }

    unsigned getMinBufferSize() const noexcept
    requires (isBuffered) {
        unsigned out = 0u;
        if constexpr (maths::BufferedExpression<TIntegrand>) out += _integrand.getMinBufferSize();
        if constexpr (maths::BufferedExpression<TJacobian>) out += _inverseJacobian.getMinBufferSize();
        return out;
    }

    void setBuffer(typename TIntegrand::Span buffer)
    requires (isBuffered) {
        std::size_t offset = 0ul;

        if constexpr (maths::BufferedExpression<TIntegrand>) {
            if constexpr (maths::BufferedExpression<TJacobian>) offset = _integrand.getMinBufferSize();
            else offset = buffer.size();
            _integrand.setBuffer({buffer.data(), offset});
        }

        if constexpr (maths::BufferedExpression<TJacobian>) {
            _inverseJacobian.setBuffer({buffer.data() + offset, buffer.data() + buffer.size()});
        }
    }

private:
    TIntegrand _integrand;

    TJacobian _inverseJacobian;
}; // class TransformedIntegrand


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
TransformedIntegrand<TIntegrand,TJacobian>
makeTransformedIntegrand(RightRef<TIntegrand> rIntegrand,
                         RightRef<TJacobian> rJacobian) {
    return TransformedIntegrand(std::move(rIntegrand), std::move(rJacobian));
}


} // namespace cie::fem
