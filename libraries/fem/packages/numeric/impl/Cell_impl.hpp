#pragma once

// help the language server
#include "packages/numeric/inc/Cell.hpp"
#include <type_traits>


namespace cie::fem {


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
CellBase<Dim,TValue,TSpatialTransform,TData>::CellBase() noexcept
    : _impl()
{
    this->ansatzID() = 0u;
    this->id() = 0ul;
}


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
CellBase<Dim,TValue,TSpatialTransform,TData>::CellBase(VertexID id,
                                                       AnsatzSpaceID ansatzID,
                                                       OrientedAxes<Dimension> axes,
                                                       RightRef<SpatialTransform> rSpatialTransform) noexcept
requires std::is_same_v<TData,void>
    : _impl(ansatzID, id, axes, std::move(rSpatialTransform), {})
{
    this->inverseSpatialTransform() = rSpatialTransform.makeInverse();
}


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
CellBase<Dim,TValue,TSpatialTransform,TData>::CellBase(VertexID id,
                                                       AnsatzSpaceID ansatzID,
                                                       OrientedAxes<Dimension> axes,
                                                       RightRef<SpatialTransform> rSpatialTransform,
                                                       typename VoidSafe<TData,int>::RightRef rData) noexcept
requires (!std::is_same_v<TData,void>)
    : _impl(ansatzID, id, axes, std::move(rSpatialTransform), {}, std::move(rData))
{
    this->inverseSpatialTransform() = rSpatialTransform.makeInverse();
}


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void CellBase<Dim,TValue,TSpatialTransform,TData>::transform(Ref<const ConstLocalSpan> in,
                                                             Ref<const GlobalSpan> out) const noexcept {
    this->spatialTransform().evaluate(
        Kernel<Dimension,Value>::decay(in),
        Kernel<Dimension,Value>::decay(out));
}


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void CellBase<Dim,TValue,TSpatialTransform,TData>::transform(Ref<const ConstGlobalSpan> in,
                                                             Ref<const LocalSpan> out) const noexcept {
    this->inverseSpatialTransform().evaluate(
        Kernel<Dimension,Value>::decay(in),
        Kernel<Dimension,Value>::decay(out));
}


template <unsigned Dim,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
typename TSpatialTransform::Derivative
CellBase<Dim,TValue,TSpatialTransform,TData>::makeJacobian() const {
    return this->spatialTransform().makeDerivative();
}


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void io::GraphML::Serializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>::header(Ref<io::GraphML::XMLElement> rElement) const {
    io::GraphML::XMLElement defaultElement = rElement.addChild("default");
    Value instance;
    this->operator()(defaultElement, instance);
}


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void io::GraphML::Serializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>::operator()(Ref<io::GraphML::XMLElement> rElement,
                                                                                             Ref<const Value> rInstance) const {
    // Serialize trivial (integral / floating point) members.
    rElement.addAttribute("id", std::to_string(rInstance.id()));
    rElement.addAttribute("iAnsatz", std::to_string(rInstance.ansatzID()));

    // Serialize custom trivial members.
    {
        std::stringstream buffer;
        buffer << rInstance.axes();
        rElement.addAttribute("axes", buffer.view());
    }

    // Serialize complex members.
    io::GraphML::Serializer<typename Value::SpatialTransform>()(
        rElement,
        rInstance.spatialTransform());

    // TData is entirely user-defined, so it must either be handled as a
    // complex member or a trivial one.
    if constexpr (std::is_integral_v<TData> || std::is_floating_point_v<TData>) {
        if constexpr (std::is_same_v<TData,bool>) {
            std::stringstream buffer;
            buffer << std::boolalpha << rInstance.data();
            rElement.addAttribute("data", buffer.view());
        } else if constexpr (std::is_integral_v<TData>) {
            rElement.addAttribute("data", std::to_string(rInstance.data()));
        } else {
            std::stringstream buffer;
            buffer << std::scientific << std::setprecision(9) << rInstance.data();
            rElement.addAttribute("data", buffer.view());
        }
    } else if constexpr (!std::is_same_v<TData,void>) {
        io::GraphML::Serializer<TData>()(
            rElement,
            rInstance.data());
    }
}


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void io::GraphML::Deserializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>::onElementBegin(Ptr<void> pThis,
                                                                                                   std::string_view elementName,
                                                                                                   std::span<io::GraphML::AttributePair> attributes) {
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    // Deserialize the cell ID.
    {
        const auto it = std::find_if(
            attributes.begin(),
            attributes.end(),
            [] (const auto pair) {
                return pair.first == "id";
            });
        char* pEnd;
        rThis.instance().id() = std::strtoul(it->second.data(), &pEnd, 10);
    }

    // Deserialize the ansatz space ID.
    {
        const auto it = std::find_if(
            attributes.begin(),
            attributes.end(),
            [] (const auto pair) {
                return pair.first == "iAnsatz";
            });
        char* pEnd;
        rThis.instance().ansatzID() = std::strtoul(it->second.data(), &pEnd, 10);
    }

    // Deserialize the cell's orientation in topological space.
    {
        const auto it = std::find_if(
            attributes.begin(),
            attributes.end(),
            [] (const auto pair) {
                return pair.first == "axes";
            });
        rThis.instance().axes() = OrientedAxes<Dimension>(it->second);
    }

    // Deserialize the spatial transformation.
    {
        using SubSerializer = io::GraphML::Deserializer<typename Value::SpatialTransform>;
        rThis.sax().push({
            SubSerializer::make(rThis.instance().spatialTransform(), rThis.sax(), elementName),
            SubSerializer::onElementBegin,
            SubSerializer::onText,
            SubSerializer::onElementEnd
        });
    }

    // Deserialize user data.
    if constexpr (std::is_integral_v<TData> || std::is_floating_point_v<TData>) {
        const auto it = std::find_if(
            attributes.begin(),
            attributes.end(),
            [] (const auto pair) {
                return pair.first == "data";
            });

        if constexpr (std::is_same_v<TData,bool>) {
            std::istringstream buffer(it->second);
            buffer >> std::boolalpha >> rThis.instance().data();
        } else if constexpr (std::is_integral_v<TData>) {
            char* pEnd;
            rThis.instance() = std::strtoul(it->second.data(), &pEnd, 10);
        } else if constexpr (std::is_floating_point_v<TData>) {
            std::istringstream buffer(it->second);
            buffer >> rThis.instance().data();
        }
    } else if constexpr (!std::is_same_v<TData,void>) {
        {
            using SubSerializer = io::GraphML::Deserializer<TData>;
            rThis.sax().push({
                SubSerializer::make(rThis.instance().data(), rThis.sax(), elementName),
                SubSerializer::onElementBegin,
                SubSerializer::onText,
                SubSerializer::onElementEnd
            });
        }
    }
}


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void io::GraphML::Deserializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>::onText(Ptr<void>,
                                                                                           std::string_view elementName) {
    CIE_THROW(
        Exception,
        "Unexpected text data while parsing a(n) " << elementName << " element in GraphML.")
}


template <unsigned Dimension,
          cie::concepts::Numeric TValue,
          maths::SpatialTransform TSpatialTransform,
          class TData>
void io::GraphML::Deserializer<CellBase<Dimension,TValue,TSpatialTransform,TData>>::onElementEnd(Ptr<void> pThis,
                                                                                                 std::string_view elementName) {
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance().inverseSpatialTransform() = rThis.instance().spatialTransform.makeInverse();
    rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem
