#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/GraphTraits.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/utilities/inc/kernel.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <tuple>


namespace cie::fem {


template <class T>
concept CellLike =
   std::is_same_v<std::remove_cvref_t<decltype(T::ParametricDimension)>,unsigned>
&& std::is_same_v<std::remove_cvref_t<decltype(T::PhysicalDimension)>,unsigned>
&& cie::concepts::Numeric<typename T::Value>
&& requires (const T& rConstInstance,
             std::span<const ParametricCoordinate<typename T::Value>,T::ParametricDimension> constParametricSpan,
             std::span<ParametricCoordinate<typename T::Value>,T::ParametricDimension> parametricSpan,
             std::span<const PhysicalCoordinate<typename T::Value>,T::PhysicalDimension> constPhysicalSpan,
             std::span<PhysicalCoordinate<typename T::Value>,T::PhysicalDimension> physicalSpan) {
    {rConstInstance.transform(constParametricSpan, physicalSpan)}   -> std::same_as<void>;      // <== transform from local to global space
    {rConstInstance.transform(constPhysicalSpan, parametricSpan)}   -> std::same_as<void>;      // <== transform from global to local space
    {rConstInstance.makeJacobian()}                                 -> maths::JacobianExpression;
    {rConstInstance.makeJacobianInverse()}                          -> maths::JacobianExpression;
    {rConstInstance.id()}                                           -> std::same_as<VertexID>;
    {rConstInstance.ansatzID()}                                     -> cie::concepts::UnsignedInteger;
    {rConstInstance.makeSpatialTransform()}                         -> maths::SpatialTransform;
}; // concept Cell


template <class T>
concept DiscretizationLike
=   GraphLike<T>
&&  CellLike<typename T::Vertex::Data>;


template <unsigned ParametricDim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData = void,
          unsigned PhysicalDim = ParametricDim>
class CellBase {
public:
    constexpr inline static unsigned ParametricDimension = ParametricDim;

    constexpr inline static unsigned PhysicalDimension = PhysicalDim;

    using ConstParametricSpan = std::span<const ParametricCoordinate<TValue>,ParametricDimension>;

    using ParametricSpan = std::span<ParametricCoordinate<TValue>,ParametricDimension>;

    using ConstPhysicalSpan = std::span<const PhysicalCoordinate<TValue>,PhysicalDimension>;

    using PhysicalSpan = std::span<PhysicalCoordinate<TValue>,PhysicalDimension>;

    using AnsatzSpaceID = unsigned short;

    using Value = TValue;

    using Data = TData;

    using SpatialTransform = TSpatialTransform;

    CellBase() noexcept;

    CellBase(
        VertexID id,
        AnsatzSpaceID ansatzID,
        OrientedAxes<ParametricDimension> axes,
        RightRef<SpatialTransform> rSpatialTransform) noexcept
    requires std::is_same_v<TData,void>;

    CellBase(
        VertexID id,
        AnsatzSpaceID ansatzID,
        OrientedAxes<ParametricDimension> axes,
        RightRef<SpatialTransform> rSpatialTransform,
        typename VoidSafe<TData,int>::RightRef rData) noexcept
    requires (!std::is_same_v<TData,void>);

    void transform(
        Ref<const ConstParametricSpan> in,
        Ref<const PhysicalSpan> out) const noexcept;

    void transform(
        Ref<const ConstPhysicalSpan> in,
        Ref<const ParametricSpan> out) const noexcept;

    typename TSpatialTransform::Derivative makeJacobian() const;

    typename TSpatialTransform::Inverse::Derivative makeJacobianInverse() const;

    [[nodiscard]] constexpr AnsatzSpaceID ansatzID() const noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr VertexID id() const noexcept {
        return std::get<1>(_impl);
    }

    [[nodiscard]] constexpr OrientedAxes<ParametricDimension> axes() const noexcept {
        return std::get<2>(_impl);
    }

    [[nodiscard]] Ref<const SpatialTransform> makeSpatialTransform() const noexcept {
        return std::get<3>(_impl);
    }

    [[nodiscard]] Ref<const typename SpatialTransform::Inverse> makeInverseSpatialTransform() const noexcept {
        return std::get<4>(_impl);
    }

    [[nodiscard]] constexpr typename VoidSafe<const TData>::Ref data() const noexcept
    requires (!std::is_same_v<TData,void>) {
        return std::get<5>(_impl);
    }

    [[nodiscard]] constexpr typename VoidSafe<TData>::Ref data() noexcept
    requires (!std::is_same_v<TData,void>) {
        return std::get<5>(_impl);
    }

protected:
    friend struct io::GraphML::Serializer<CellBase>;

    friend struct io::GraphML::Deserializer<CellBase>;

    [[nodiscard]] constexpr Ref<AnsatzSpaceID> ansatzID() noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr Ref<VertexID> id() noexcept {
        return std::get<1>(_impl);
    }

    [[nodiscard]] constexpr Ref<OrientedAxes<ParametricDimension>> axes() noexcept {
        return std::get<2>(_impl);
    }

    [[nodiscard]] constexpr Ref<const TSpatialTransform> spatialTransform() const noexcept {
        return std::get<3>(_impl);
    }

    [[nodiscard]] constexpr Ref<TSpatialTransform> spatialTransform() noexcept {
        return std::get<3>(_impl);
    }

    [[nodiscard]] constexpr Ref<const typename TSpatialTransform::Inverse> inverseSpatialTransform() const noexcept {
        return std::get<4>(_impl);
    }

    [[nodiscard]] constexpr Ref<typename TSpatialTransform::Inverse> inverseSpatialTransform() noexcept {
        return std::get<4>(_impl);
    }

    using Impl = std::conditional_t<
        std::is_same_v<TData,void>,
        std::tuple<
            AnsatzSpaceID,
            VertexID,
            OrientedAxes<ParametricDimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse>,
        std::tuple<
            AnsatzSpaceID,
            VertexID,
            OrientedAxes<ParametricDimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse,
            TData>
    >;
    Impl _impl;
}; // class CellBase


template <
    unsigned ParametricDimension,
    cie::concepts::Numeric TValue,
    maths::SpatialTransform TSpatialTransform,
    class TData,
    unsigned PhysicalDimension>
struct io::GraphML::Serializer<CellBase<ParametricDimension,TValue,TSpatialTransform,TData,PhysicalDimension>> {
    using Value = CellBase<ParametricDimension,TValue,TSpatialTransform,TData,PhysicalDimension>;

    void header(Ref<io::GraphML::XMLElement> rElement) const;

    void operator()(
        Ref<io::GraphML::XMLElement> rElement,
        Ref<const Value> rInstance) const;
}; // struct io::GraphML::Serializer<CellBase>


template <
    unsigned ParametricDimension,
    cie::concepts::Numeric TValue,
    maths::SpatialTransform TSpatialTransform,
    class TData,
    unsigned PhysicalDimension>
struct io::GraphML::Deserializer<CellBase<ParametricDimension,TValue,TSpatialTransform,TData,PhysicalDimension>>
    : public io::GraphML::DeserializerBase<CellBase<ParametricDimension,TValue,TSpatialTransform,TData,PhysicalDimension>>
{
    using Value = CellBase<ParametricDimension,TValue,TSpatialTransform,TData>;

    using io::GraphML::DeserializerBase<Value>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair> attributes);

    static void onText(Ptr<void>, std::string_view);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);
}; // struct io::GraphML::Deserializer<CellBase>


} // namespace cie::fem

#include "packages/numeric/impl/Cell_impl.hpp"
