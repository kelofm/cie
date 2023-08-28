#ifndef CIE_LINALG_QUATERNION_IMPL_HPP
#define CIE_LINALG_QUATERNION_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/overloads/inc/vectoroperators.hpp"
#include "packages/types/inc/typeoperations.hpp"


namespace cie::linalg {


template <concepts::Numeric NT>
Quaternion<NT>::Quaternion(typename Quaternion<NT>::component_container&& r_components) :
    _components(std::move(r_components))
{
}


template <concepts::Numeric NT>
Quaternion<NT>::Quaternion(const typename Quaternion<NT>::component_container& r_components) :
    _components(r_components)
{
}


template <concepts::Numeric NT>
Quaternion<NT>::Quaternion() :
    _components({1, 0, 0, 0})
{
}


template <concepts::Numeric NT>
inline void
Quaternion<NT>::conjugate()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it_component = this->_components.begin() + 1;
    *it_component++ *= -1;
    *it_component++ *= -1;
    *it_component *= -1;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
void
Quaternion<NT>::normalize()
{
    CIE_BEGIN_EXCEPTION_TRACING

    NT magnitude = std::sqrt(this->normSquared());
    this->_components /= magnitude;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
inline NT
Quaternion<NT>::normSquared() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it_component = this->_components.begin();
    NT output = (*it_component) * (*it_component);
    ++it_component;
    output += (*it_component) * (*it_component);
    ++it_component;
    output += (*it_component) * (*it_component);
    ++it_component;
    output += (*it_component) * (*it_component);

    return output;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
template <class TAxis>
requires concepts::Container<TAxis,NT>
void
Quaternion<NT>::loadFromAxisAndAngle(Ref<const TAxis> r_axis, NT angle)
{
    CIE_BEGIN_EXCEPTION_TRACING

    NT halfAngle   = angle / 2;
    NT coefficient = std::sin(halfAngle);

    NT axisNorm = norm(r_axis);

    this->_components[0] = std::cos(halfAngle);
    this->_components[1] = coefficient * r_axis[0] / axisNorm;
    this->_components[2] = coefficient * r_axis[1] / axisNorm;
    this->_components[3] = coefficient * r_axis[2] / axisNorm;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
template <class MatrixType>
void
Quaternion<NT>::toRotationMatrix(MatrixType& r_matrix) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK(3 <= r_matrix.numberOfRows() && 3 <= r_matrix.numberOfColumns())

    NT normSquared = this->normSquared();
    CIE_DIVISION_BY_ZERO_CHECK(normSquared != 0)

    const NT w = this->_components[0];
    const NT i = this->_components[1];
    const NT j = this->_components[2];
    const NT k = this->_components[3];

    const NT ii = i * i / normSquared;
    const NT jj = j * j / normSquared;
    const NT kk = k * k / normSquared;

    const NT ij = i * j / normSquared;
    const NT ik = i * k / normSquared;
    const NT jk = j * k / normSquared;

    const NT wi = w * i / normSquared;
    const NT wj = w * j / normSquared;
    const NT wk = w * k / normSquared;

    r_matrix(0, 0) = 1 - 2*(jj + kk);
    r_matrix(0, 1) = 2*(ij - wk);
    r_matrix(0, 2) = 2*(ik + wj);

    r_matrix(1, 0) = 2*(ij + wk);
    r_matrix(1, 1) = 1 - 2*(ii + kk);
    r_matrix(1, 2) = 2*(jk - wi);

    r_matrix(2, 0) = 2*(ik - wj);
    r_matrix(2, 1) = 2*(jk + wi);
    r_matrix(2, 2) = 1 - 2*(ii + jj);

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
inline const typename Quaternion<NT>::component_container&
Quaternion<NT>::components() const
{
    return this->_components;
}


template <concepts::Numeric NT>
inline NT
Quaternion<NT>::operator[](Size index) const
{
    CIE_OUT_OF_RANGE_CHECK(index < this->_components.size())
    return _components[index];
}


template <concepts::Numeric NT>
inline NT&
Quaternion<NT>::operator[](Size index)
{
    CIE_OUT_OF_RANGE_CHECK(index < this->_components.size())
    return _components[index];
}


template <concepts::Numeric NT>
inline typename Quaternion<NT>::const_iterator
Quaternion<NT>::begin() const
{
    return _components.begin();
}


template <concepts::Numeric NT>
inline typename Quaternion<NT>::iterator
Quaternion<NT>::begin()
{
    return _components.begin();
}


template <concepts::Numeric NT>
inline typename Quaternion<NT>::const_iterator
Quaternion<NT>::end() const
{
    return _components.end();
}


template <concepts::Numeric NT>
inline typename Quaternion<NT>::iterator
Quaternion<NT>::end()
{
    return _components.end();
}


template <concepts::Numeric NT>
inline void
Quaternion<NT>::operator*=(NT coefficient)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using ::operator*=;
    this->_components *= coefficient;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
inline void
Quaternion<NT>::operator/=(NT denominator)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using ::operator/=;
    this->_components /= denominator;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
inline void
Quaternion<NT>::operator+=(const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using ::operator+=;
    this->_components += r_rhs._components;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
inline void
Quaternion<NT>::operator-=(const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using ::operator-=;
    this->_components -= r_rhs._components;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
void
Quaternion<NT>::operator*=(const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto tmp = this->_components;

    this->_components[0] = tmp[0]*r_rhs[0] - tmp[1]*r_rhs[1] - tmp[2]*r_rhs[2] - tmp[3]*r_rhs[3];
    this->_components[1] = tmp[0]*r_rhs[1] + tmp[1]*r_rhs[0] + tmp[2]*r_rhs[3] - tmp[3]*r_rhs[2];
    this->_components[2] = tmp[0]*r_rhs[2] + tmp[2]*r_rhs[0] + tmp[3]*r_rhs[1] - tmp[1]*r_rhs[3];
    this->_components[3] = tmp[0]*r_rhs[3] + tmp[3]*r_rhs[0] + tmp[1]*r_rhs[2] - tmp[2]*r_rhs[1];

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
Quaternion<NT> operator+(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    Quaternion<NT> output(r_lhs);
    output += r_rhs;
    return output;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Numeric NT>
Quaternion<NT> operator-(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    Quaternion<NT> output(r_lhs);
    output -= r_rhs;
    return output;

    CIE_END_EXCEPTION_TRACING
}

template <concepts::Numeric NT>
Quaternion<NT> operator*(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    Quaternion<NT> output(r_lhs);
    output *= r_rhs;
    return output;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::linalg

#endif