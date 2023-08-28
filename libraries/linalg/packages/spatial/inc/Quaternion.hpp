#ifndef LINALG_QUATERNION_HPP
#define LINALG_QUATERNION_HPP

// ---- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <array>


namespace cie::linalg {


template <concepts::Numeric NT>
class Quaternion
{
public:
    using component_container = StaticArray<NT,4>;
    using iterator            = component_container::iterator;
    using const_iterator      = component_container::const_iterator;

public:
    Quaternion(component_container&& r_components);

    Quaternion(const component_container& r_components);

    Quaternion();

    Quaternion(Quaternion<NT>&& r_rhs) = default;

    Quaternion(const Quaternion<NT>& r_rhs) = default;

    Quaternion<NT>& operator=(const Quaternion<NT>& r_rhs) = default;

public:
    void conjugate();

    void normalize();

    NT normSquared() const;

    template <class TAxis>
    requires concepts::Container<TAxis,NT>
    void loadFromAxisAndAngle(Ref<const TAxis> r_axis, NT angle);

    template <class MatrixType>
    void toRotationMatrix(MatrixType& r_matrix) const;

    const component_container& components() const;

    NT operator[](Size index) const;
    NT& operator[](Size index);

    const_iterator begin() const;
    iterator begin();

    const_iterator end() const;
    iterator end();

public:
    void operator*=(NT coefficient);

    void operator/=(NT denominator);

    void operator+=(const Quaternion<NT>& r_rhs);

    void operator-=(const Quaternion<NT>& r_rhs);

    void operator*=(const Quaternion<NT>& r_rhs);

protected:
    component_container  _components;
};


template <concepts::Numeric NT>
Quaternion<NT> operator+(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs);

template <concepts::Numeric NT>
Quaternion<NT> operator-(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs);

template <concepts::Numeric NT>
Quaternion<NT> operator*(const Quaternion<NT>& r_lhs, const Quaternion<NT>& r_rhs);


} // namespace cie::linalg

#include "packages/spatial/impl/Quaternion_impl.hpp"

#endif