#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
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
   std::is_same_v<std::remove_cvref_t<decltype(T::Dimension)>,unsigned>
&& cie::concepts::Numeric<typename T::Value>
&& requires (const T& rConstInstance,
             std::span<const typename Kernel<T::Dimension,typename T::Value>::LocalCoordinate,T::Dimension> constLocalSpan,
             std::span<typename Kernel<T::Dimension,typename T::Value>::LocalCoordinate,T::Dimension> localSpan,
             std::span<const typename Kernel<T::Dimension,typename T::Value>::GlobalCoordinate,T::Dimension> constGlobalSpan,
             std::span<typename Kernel<T::Dimension,typename T::Value>::GlobalCoordinate,T::Dimension> globalSpan) {
    {rConstInstance.transform(constLocalSpan, globalSpan)}  -> std::same_as<void>;      // <== transform from local to global space
    {rConstInstance.transform(constGlobalSpan, localSpan)}  -> std::same_as<void>;      // <== transform from global to local space
    {rConstInstance.makeJacobian()}                         -> maths::JacobianExpression;
    {rConstInstance.id()}                                   -> std::same_as<VertexID>;
    {rConstInstance.ansatzSpaceID()}                        -> ::cie::concepts::UnsignedInteger;
}; // concept Cell


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData = void>
class CellBase {
public:
    constexpr inline static unsigned Dimension = Dim;

    using ConstLocalSpan = std::span<const typename Kernel<Dimension,TValue>::LocalCoordinate,Dimension>;

    using LocalSpan = std::span<typename Kernel<Dimension,TValue>::LocalCoordinate,Dimension>;

    using ConstGlobalSpan = std::span<const typename Kernel<Dimension,TValue>::GlobalCoordinate,Dimension>;

    using GlobalSpan = std::span<typename Kernel<Dimension,TValue>::GlobalCoordinate,Dimension>;

    using AnsatzSpaceID = unsigned short;

    using Value = TValue;

    using Data = TData;

    using SpatialTransform = TSpatialTransform;

    CellBase() noexcept;

    CellBase(VertexID id,
             AnsatzSpaceID ansatzSpaceID,
             OrientedAxes<Dimension> axes,
             RightRef<SpatialTransform> rSpatialTransform) noexcept
    requires std::is_same_v<TData,void>;

    CellBase(VertexID id,
             AnsatzSpaceID ansatzSpaceID,
             OrientedAxes<Dimension> axes,
             RightRef<SpatialTransform> rSpatialTransform,
             typename VoidSafe<TData,int>::RightRef rData) noexcept
    requires (!std::is_same_v<TData,void>);

    void transform(Ref<const ConstLocalSpan> in, Ref<const GlobalSpan> out) const noexcept;

    void transform(Ref<const ConstGlobalSpan> in, Ref<const LocalSpan> out) const noexcept;

    typename SpatialTransform::Derivative makeJacobian() const;

    [[nodiscard]] constexpr AnsatzSpaceID ansatzSpaceID() const noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr VertexID id() const noexcept {
        return std::get<1>(_impl);
    }

    [[nodiscard]] constexpr OrientedAxes<Dimension> axes() const noexcept {
        return std::get<2>(_impl);
    }

    [[nodiscard]] constexpr typename VoidSafe<const TData>::Ref data() const noexcept
    requires (!std::is_same_v<TData,void>) {
        return std::get<5>(_impl);
    }

    [[nodiscard]] constexpr typename VoidSafe<TData>::Ref data() noexcept
    requires (!std::is_same_v<TData,void>) {
        return std::get<5>(_impl);
    }

private:
    friend struct io::GraphML::Serializer<CellBase>;

    friend struct io::GraphML::Deserializer<CellBase>;

    [[nodiscard]] constexpr Ref<AnsatzSpaceID> ansatzSpaceID() noexcept {
        return std::get<0>(_impl);
    }

    [[nodiscard]] constexpr Ref<VertexID> id() noexcept {
        return std::get<1>(_impl);
    }

    [[nodiscard]] constexpr Ref<OrientedAxes<Dimension>> axes() noexcept {
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
            OrientedAxes<Dimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse>,
        std::tuple<
            AnsatzSpaceID,
            VertexID,
            OrientedAxes<Dimension>,
            TSpatialTransform,
            typename TSpatialTransform::Inverse,
            TData>
    >;
    Impl _impl;
}; // class CellBase


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
struct io::GraphML::Serializer<CellBase<Dimension,TValue,TSpatialTransform,TData>> {
    using Value = CellBase<Dimension,TValue,TSpatialTransform,TData>;

    void header(Ref<io::GraphML::XMLElement> rElement) const;

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const Value> rInstance) const;
}; // struct io::GraphML::Serializer<CellBase>


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
struct io::GraphML::Deserializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>
    : public io::GraphML::DeserializerBase<CellBase<Dimension,TValue,TSpatialTransform,TData>>
{
    using Value = CellBase<Dimension,TValue,TSpatialTransform,TData>;

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
