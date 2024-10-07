#ifndef CIE_FEM_LINEAR_ISOTROPIC_STIFFNESS_INTEGRAND_HPP
#define CIE_FEM_LINEAR_ISOTROPIC_STIFFNESS_INTEGRAND_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <memory> // shared_ptr
#include <span> // span


namespace cie::fem::maths {


template <Expression TAnsatzDerivatives>
class LinearIsotropicStiffnessIntegrand : public ExpressionTraits<typename TAnsatzDerivatives::Value>
{
public:
    static constexpr unsigned Dimension = TAnsatzDerivatives::Dimension;

    using typename ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename ExpressionTraits<Value>::ConstIterator;

    using typename ExpressionTraits<Value>::Iterator;

public:
    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<const std::shared_ptr<const TAnsatzDerivatives>> pAnsatzDerivatives);

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<const std::shared_ptr<const TAnsatzDerivatives>> pAnsatzDerivatives,
                                      std::span<Value> buffer);

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _modulus;

    std::shared_ptr<const TAnsatzDerivatives> _pAnsatzDerivatives;

    std::span<Value> _buffer;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem::maths

#include "packages/maths/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"

#endif
