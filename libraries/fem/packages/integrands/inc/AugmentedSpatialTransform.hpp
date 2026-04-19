#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/GraphML.hpp"


namespace cie::fem {


template <maths::SpatialTransform T>
class AugmentedSpatialTransform : public T {
public:
    using typename T::Span;

    using typename T::ConstSpan;

    using typename T::BufferSpan;

    using T::T;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const {
            ConstSpan restrictedIn(
                in.data(),
                T::ParametricDimension);
            Span restrictedOut(
                out.data(),
                T::PhysicalDimension);
            T::evaluate(
                restrictedIn,
                restrictedOut,
                buffer);
            std::copy(
                in.begin() + T::ParametricDimension,
                in.end(),
                out.begin() + T::PhysicalDimension);
    }
}; // class AugmentedSpatialTransform


} // namespace cie::fem


namespace cie::fem::io {


template <class T>
struct GraphML::Serializer<AugmentedSpatialTransform<T>>
    : public GraphML::Serializer<T>
{};


template <class T>
struct GraphML::Deserializer<AugmentedSpatialTransform<T>>
    : public GraphML::Deserializer<T> {
        using GraphML::Deserializer<T>::Deserializer;
}; // GraphML::Deserializer<ScaleTranslateTransform>


} // namespace cie::fem::io
