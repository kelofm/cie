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
#include "packages/stl_extension/inc/StrongTypeDef.hpp"

// --- STL Includes ---
#include <tuple>


namespace cie::fem {


/// @brief Ansatz space identifier type.
/// @class AnsatzID
/// @see StrongTypeDef
CIE_STRONG_TYPEDEF(unsigned, AnsatzID);


template <class T>
concept CellLike
=  std::is_same_v<std::remove_cvref_t<decltype(T::ParametricDimension)>,unsigned>
&& std::is_same_v<std::remove_cvref_t<decltype(T::PhysicalDimension)>,unsigned>
&& cie::concepts::Numeric<typename T::Value>
&& std::is_same_v<typename T::BufferSpan,typename maths::ExpressionTraits<typename T::Value>::BufferSpan>
&& requires (const T& rConstInstance,
             std::span<const ParametricCoordinate<typename T::Value>,T::ParametricDimension> constParametricSpan,
             std::span<ParametricCoordinate<typename T::Value>,T::ParametricDimension> parametricSpan,
             std::span<const PhysicalCoordinate<typename T::Value>,T::PhysicalDimension> constPhysicalSpan,
             std::span<PhysicalCoordinate<typename T::Value>,T::PhysicalDimension> physicalSpan,
             typename T::BufferSpan bufferSpan) {
    {rConstInstance.transform(constParametricSpan, physicalSpan, bufferSpan)}   -> std::same_as<void>;      // <== transform from local to global space
    {rConstInstance.transform(constPhysicalSpan, parametricSpan, bufferSpan)}   -> std::same_as<void>;      // <== transform from global to local space
    {rConstInstance.makeJacobian()}                                             -> maths::JacobianExpression;
    {rConstInstance.makeJacobianInverse()}                                      -> maths::JacobianExpression;
    {rConstInstance.id()}                                                       -> std::same_as<VertexID>;
    {rConstInstance.ansatzID()}                                                 -> std::same_as<AnsatzID>;
    {rConstInstance.makeSpatialTransform()}                                     -> maths::SpatialTransform;
}; // concept CellLike


template <class T>
concept CellBoundaryLike
=  std::is_same_v<std::remove_cvref_t<decltype(T::Dimension)>,unsigned>
&& requires (const T& rConstInstance) {
    {rConstInstance.boundary()} -> std::same_as<BoundaryID>;
}; // concept CellBoundaryLike


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

    using Value = TValue;

    using Data = TData;

    using SpatialTransform = TSpatialTransform;

    using BufferSpan = typename maths::ExpressionTraits<TValue>::BufferSpan;

    CellBase() noexcept;

    CellBase(
        VertexID id,
        AnsatzID ansatzID,
        OrientedAxes<ParametricDimension> axes,
        RightRef<SpatialTransform> rSpatialTransform) noexcept
    requires std::is_same_v<TData,void>;

    CellBase(
        VertexID id,
        AnsatzID ansatzID,
        OrientedAxes<ParametricDimension> axes,
        RightRef<SpatialTransform> rSpatialTransform,
        typename VoidSafe<TData,int>::RightRef rData) noexcept
    requires (!std::is_same_v<TData,void>);

    void transform(
        Ref<const ConstParametricSpan> in,
        Ref<const PhysicalSpan> out,
        Ref<const BufferSpan> buffer) const noexcept;

    void transform(
        Ref<const ConstPhysicalSpan> in,
        Ref<const ParametricSpan> out,
        Ref<const BufferSpan> buffer) const noexcept;

    [[nodiscard]] typename TSpatialTransform::Derivative makeJacobian() const;

    [[nodiscard]] typename TSpatialTransform::Inverse::Derivative makeJacobianInverse() const;

    [[nodiscard]] constexpr VertexID id() const noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr AnsatzID ansatzID() const noexcept {
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

    [[nodiscard]] constexpr Ref<VertexID> id() noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr Ref<AnsatzID> ansatzID() noexcept {
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
            VertexID,
            AnsatzID,
            OrientedAxes<ParametricDimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse>,
        std::tuple<
            VertexID,
            AnsatzID,
            OrientedAxes<ParametricDimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse,
            TData>
    >;
    Impl _impl;
}; // class CellBase


template <CellLike TCell>
class IndirectCell {
public:
    constexpr inline static unsigned ParametricDimension = TCell::ParametricDimension;

    constexpr inline static unsigned PhysicalDimension = TCell::PhysicalDimension;

    using Value = typename TCell::Value;

    using Data = typename TCell::Data;

    using ConstParametricSpan = std::span<const ParametricCoordinate<Value>,ParametricDimension>;

    using ParametricSpan = std::span<ParametricCoordinate<Value>,ParametricDimension>;

    using ConstPhysicalSpan = std::span<const PhysicalCoordinate<Value>,PhysicalDimension>;

    using PhysicalSpan = std::span<PhysicalCoordinate<Value>,PhysicalDimension>;

    using SpatialTransform = typename TCell::SpatialTransform;

    using BufferSpan = typename TCell::BufferSpan;

    constexpr IndirectCell() noexcept
        : _pCell(nullptr)
    {}

    constexpr IndirectCell(Ref<TCell> rCell) noexcept
        : _pCell(&rCell)
    {}

    void transform(
        Ref<const ConstParametricSpan> in,
        Ref<const PhysicalSpan> out,
        Ref<const BufferSpan> buffer) const noexcept {
            _pCell->transform(in, out, buffer);
    }

    void transform(
        Ref<const ConstPhysicalSpan> in,
        Ref<const ParametricSpan> out,
        Ref<const BufferSpan> buffer) const noexcept {
            _pCell->transform(in, out, buffer);
    }

    [[nodiscard]] typename SpatialTransform::Derivative makeJacobian() const {
        return _pCell->makeJacobian();
    }

    [[nodiscard]] typename SpatialTransform::Inverse::Derivative makeJacobianInverse() const {
        return _pCell->makeJacobianInverse();
    }

    [[nodiscard]] constexpr VertexID id() const noexcept {
        return _pCell->id();
    }

    [[nodiscard]] constexpr OrientedAxes<ParametricDimension> axes() const noexcept {
        return _pCell->axes();
    }

    [[nodiscard]] Ref<const SpatialTransform> makeSpatialTransform() const noexcept {
        return _pCell->makeSpatialTransform();
    }

    [[nodiscard]] Ref<const typename SpatialTransform::Inverse> makeInverseSpatialTransform() const noexcept {
        return _pCell->makeInverseSpatialTransform();
    }

    [[nodiscard]] constexpr typename VoidSafe<const Data>::Ref data() const noexcept
    requires (!std::is_same_v<Data,void>) {
        return _pCell->data();
    }

    [[nodiscard]] constexpr typename VoidSafe<Data>::Ref data() noexcept
    requires (!std::is_same_v<Data,void> && !std::is_const_v<TCell>) {
        return _pCell->data();
    }

private:
    Ptr<TCell> _pCell;
}; // class IndirectCell


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

#include "packages/numeric/impl/CellBase_impl.hpp"
