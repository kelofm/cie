#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/BoxBoundable.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"


namespace cie::fem {


/// @brief Data structure unique to each @ref Graph::Vertex "cell".
struct CellData : public geo::BoxBoundable<Dimension,Scalar> {
    using geo::BoxBoundable<Dimension,Scalar>::BoundingBox;

    using LocalCoordinate = Kernel<Dimension,Scalar>::LocalCoordinate;

    using GlobalCoordinate = Kernel<Dimension,Scalar>::GlobalCoordinate;

    VertexID id;

    /// @brief Index of the cell's ansatz space in @ref MeshData::ansatzSpaces.
    unsigned short iAnsatz;

    /// @brief Local diffusivity coefficient.
    /// @details Assumed to be constant throughout the cell for simplicity.
    Scalar diffusivity;

    /// @brief Local axis orientation of the cell to enforce continuity.
    /// @details Adjacent cells must produce identical state fields
    ///          on both sides of the shared boundary. A simple way of
    ///          ensuring this is to orient the local axes of each cell
    ///          such that matching axes point in the same direction on
    ///          both sides, except the shared boundary's normal axes,
    ///          that are pointing toward each other.
    ///          @code
    ///          + ----------- +   + ----------- +
    ///          |      y      |   |      y      |
    ///          |      ^      |   |      ^      |
    ///          |      |      |   |      |      |
    ///          |      + -> x |   | x <- +      |
    ///          |             |   |             |
    ///          |             |   |             |
    ///          + ----------- +   + ----------- +
    ///          @endcode
    ///          The oriented axes define an intermediate space between
    ///          local and global space (referred to as topological space).
    ///          DoFs that don't vanish on the shared boundaries can be
    ///          merged on both sides.
    OrientedAxes<Dimension> axes;

    /// @brief Spatial transform from local to global space.
    SpatialTransform spatialTransform;

    SpatialTransform::Inverse inverseSpatialTransform;

    CellData() noexcept = default;

    CellData(VertexID id_,
             unsigned short iAnsatz_,
             Scalar diffusivity_,
             OrientedAxes<Dimension> axes_,
             RightRef<SpatialTransform> rSpatialTransform) noexcept
        : id(id_),
          iAnsatz(iAnsatz_),
          diffusivity(diffusivity_),
          axes(axes_),
          spatialTransform(std::move(rSpatialTransform)) {
        this->inverseSpatialTransform = spatialTransform.makeInverse();
    }

    void transform(Ref<const std::span<const LocalCoordinate,Dimension>> rLocalCoordinates,
                   Ref<const std::span<GlobalCoordinate,Dimension>> rGlobalCoordinates) const noexcept {
        spatialTransform.evaluate(
            Kernel<Dimension,Scalar>::decay(rLocalCoordinates),
            Kernel<Dimension,Scalar>::decay(rGlobalCoordinates));
    }

    void transform(Ref<const std::span<const GlobalCoordinate,Dimension>> rGlobalCoordinates,
                   Ref<const std::span<LocalCoordinate,Dimension>> rLocalCoordinates) const noexcept {
        spatialTransform.evaluate(
            Kernel<Dimension,Scalar>::decay(rGlobalCoordinates),
            Kernel<Dimension,Scalar>::decay(rLocalCoordinates));
    }

    constexpr VertexID ID() const noexcept {
        return id;
    }

    constexpr unsigned ansatzSpaceID() const noexcept {
        return iAnsatz;
    }

    bool at(geo::BoxBoundable<Dimension,Scalar>::Point point) const {
        const utils::Comparison<Scalar> comparison(1e-8, 1e-10);
        StaticArray<Scalar,Dimension> local;

        this->inverseSpatialTransform.evaluate(point, local);

        return std::all_of(
            local.begin(),
            local.end(),
            [&comparison](Scalar coordinate) {
                return comparison.less(std::abs(coordinate), static_cast<Scalar>(1))
                    || comparison.equal(std::abs(coordinate), static_cast<Scalar>(1));}
        );
    }

protected:
    void computeBoundingBoxImpl(BoundingBox& rBox) noexcept override
    {
        BoundingBox::Point opposite;
        StaticArray<LocalCoordinate,Dimension> localCorner;
        StaticArray<GlobalCoordinate,Dimension> globalCorner;

        std::fill_n(
            rBox.base().data(),
            Dimension,
            std::numeric_limits<Scalar>::max());
        std::fill_n(
            opposite.data(),
            Dimension,
            std::numeric_limits<Scalar>::lowest());

        StaticArray<std::uint8_t,Dimension> state {0u, 0u};
        StaticArray<LocalCoordinate,2> ordinates {-1.0, 1.0};

        do {
            // Compute the corner in local space.
            std::transform(state.begin(),
                           state.end(),
                           localCorner.begin(),
                           [&ordinates](std::uint8_t iOrdinate) {
                                return ordinates[iOrdinate];
                           });

            // Transform the corner to global space.
            this->transform(
                Kernel<Dimension,Scalar>::view(localCorner),
                Kernel<Dimension,Scalar>::view(globalCorner));

            // Extend box definition.
            for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                rBox.base()[iDimension] = std::min<Scalar>(rBox.base()[iDimension], globalCorner[iDimension]);
                opposite[iDimension] = std::max<Scalar>(opposite[iDimension], globalCorner[iDimension]);
            } // for iDimension in range(Dimension)
        } while (cie::maths::OuterProduct<Dimension>::next(2u, state.data()));

        std::transform(opposite.begin(),
                       opposite.end(),
                       rBox.base().begin(),
                       rBox.lengths().begin(),
                       std::minus<Scalar>());
    }
}; // struct CellData


template <>
struct io::GraphML::Serializer<CellData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const {
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        CellData instance;
        this->operator()(defaultElement, instance);
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const CellData> rInstance) const {
        rElement.addAttribute("iAnsatz", std::to_string(rInstance.iAnsatz));
        rElement.addAttribute("diffusivity", std::to_string(rInstance.diffusivity));

        std::stringstream buffer;
        buffer << rInstance.axes;
        rElement.addAttribute("axes", buffer.view());

        io::GraphML::Serializer<SpatialTransform>()(rElement, rInstance.spatialTransform);
    }
}; // struct GraphML::Serializer<CellData>


template <>
struct io::GraphML::Deserializer<CellData>
    : public io::GraphML::DeserializerBase<CellData>
{
    using io::GraphML::DeserializerBase<CellData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair> attributes) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "iAnsatz";
                                         });
            char* pEnd;
            rThis.instance().iAnsatz = std::strtoul(it->second.data(), &pEnd, 10);
        }

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "diffusivity";
                                         });
            char* pEnd;
            rThis.instance().diffusivity = std::strtod(it->second.data(), &pEnd);
        }

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "axes";
                                         });
            rThis.instance().axes = OrientedAxes<Dimension>(it->second);
        }

        {
            using SubSerializer = io::GraphML::Deserializer<SpatialTransform>;
            rThis.sax().push({
                SubSerializer::make(rThis.instance().spatialTransform, rThis.sax(), elementName),
                SubSerializer::onElementBegin,
                SubSerializer::onText,
                SubSerializer::onElementEnd
            });
        }

        CIE_THROW(NotImplementedException, "")
    }

    static void onText(Ptr<void>, std::string_view) noexcept {}

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName) {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.instance().inverseSpatialTransform = rThis.instance().spatialTransform.makeInverse();
        rThis.template release<Deserializer>(&rThis, elementName);
    }

}; // struct GraphML::Deserializer<CellData>


} // namespace cie::fem
