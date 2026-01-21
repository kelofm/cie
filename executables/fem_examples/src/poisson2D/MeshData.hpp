#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"
#include "poisson2D/CellData.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/QuadraturePointFactory.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/sycl.hpp"


namespace cie::fem {


class SYCLSingleton {
public:
    #ifdef CIE_ENABLE_SYCL
        static Ref<sycl::device> getDevice() {
            SYCLSingleton::init();
            return *_maybeImpl.value().pDevice;
        }

        static Ref<sycl::queue> getQueue() {
            SYCLSingleton::init();
            return *_maybeImpl.value().pQueue;
        }

        template <class T>
        static sycl::usm_allocator<T,sycl::usm::alloc::shared> makeSharedAllocator() {
            return sycl::usm_allocator<T,sycl::usm::alloc::shared>(SYCLSingleton::getQueue());
        }
    #else
        template <class T>
        static std::allocator<T> makeSharedAllocator() {
            return std::allocator<T>();
        }
    #endif

private:
    #ifdef CIE_ENABLE_SYCL
        static void init() {
            if (!_maybeImpl.has_value()) {
                auto pDevice = std::make_unique<sycl::device>(sycl::default_selector_v);
                auto pQueue = std::make_unique<sycl::queue>(*pDevice);
                _maybeImpl.emplace(std::move(pDevice), std::move(pQueue));
            }
        }

        struct Impl {
            ~Impl() {
                pQueue.reset();
                pDevice.reset();
            }
            std::unique_ptr<sycl::device> pDevice;
            std::unique_ptr<sycl::queue> pQueue;
        };

        static std::optional<Impl> _maybeImpl;
    #endif
}; // class SYCLSingleton


#ifdef CIE_ENABLE_SYCL
    std::optional<SYCLSingleton::Impl> SYCLSingleton::_maybeImpl = {};
#endif


/// @brief Data structure common to the entire @ref Graph "mesh".
class MeshData {
public:
    MeshData()
        : _ansatzSpace(),
          _ansatzDerivative(),
          _quadraturePointSet(),
          _buffer()
    {}

    MeshData(RightRef<Ansatz> rAnsatzSpace,
             RightRef<Ansatz::Derivative> rAnsatzDerivative)
        : _ansatzSpace(std::move(rAnsatzSpace)),
          _ansatzDerivative(std::move(rAnsatzDerivative)),
          _quadraturePointSet()
    {
        // Cache quadrature points.
        CellData dummyElement;

        // Generate 1D quadrature points.
        DynamicArray<QuadraturePoint<1,Scalar>> basePoints;
        OuterProductQuadraturePointFactory<Dimension,Scalar> generator;

        {
            GaussLegendreQuadrature<Scalar> quadrature(integrationOrder);
            basePoints.reserve(quadrature.numberOfNodes());
            for (std::size_t iNode=0ul; iNode<quadrature.numberOfNodes(); ++iNode) {
                basePoints.emplace_back(
                    quadrature.nodes()[iNode],
                    quadrature.weights()[iNode]);
            }
            generator = OuterProductQuadraturePointFactory<Dimension,Scalar>(basePoints);
        }

        // Generate nD quadrature points.
        constexpr std::size_t initialSize = 0x10;
        constexpr double growFactor = 2.0;
        std::size_t pointCount = 0ul;
        _quadraturePointSet.resize(initialSize);

        while (true) {
            // Extend the container if necessary.
            if (!(pointCount < _quadraturePointSet.size())) {
                _quadraturePointSet.resize(growFactor * _quadraturePointSet.size());
            }

            std::span<QuadraturePoint<Dimension,Scalar>> targetSpan(
                _quadraturePointSet.data() + pointCount,
                _quadraturePointSet.data() + _quadraturePointSet.size());

            // Request a new batch of quadrature points.
            const auto newPointCount = generator.generate(dummyElement, targetSpan);
            pointCount += newPointCount;

            if (!newPointCount) break;
        }

        _quadraturePointSet.resize(pointCount);
    }

    Ref<const Ansatz> ansatzSpace() const noexcept {
        return _ansatzSpace;
    }

    Ref<const Ansatz::Derivative> ansatzDerivative() const noexcept {
        return _ansatzDerivative;
    }

    CachedQuadraturePointFactory<Dimension,Scalar> makeQuadratureRule() const {
        return CachedQuadraturePointFactory<Dimension,Scalar>(
            std::span<const QuadraturePoint<Dimension,Scalar>>(
                _quadraturePointSet.data(),
                _quadraturePointSet.size())
        );
    }

private:
    friend struct io::GraphML::Serializer<MeshData>;

    friend struct io::GraphML::Deserializer<MeshData>;

    /// @brief Collection of all ansatz spaces the contained cells can refer to.
    Ansatz _ansatzSpace;

    /// @brief Collection of all ansatz spaces' derivatives the contained cells can refer to.
    Ansatz::Derivative _ansatzDerivative;

    /// @brief Set of quadrature points for a default local hypercube.
    /// @details These quadrature points are used while constructing
    ///          cell-specific quadrature rules.
    DynamicArray<QuadraturePoint<Dimension,Scalar>> _quadraturePointSet;

    DynamicArray<Scalar> _buffer;
}; // class MeshData


/// @brief Serializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Serializer<MeshData> {
    void header(Ref<io::GraphML::XMLElement> rElement) const {
        // Add default value to the header.
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        MeshData instance;
        this->operator()(defaultElement, instance);

        // Add description to the header.
        io::GraphML::XMLElement descriptionElement = rElement.addChild("desc");
        std::stringstream description;
        description << "Data structure shared by all cells and boundaries of the mesh. "
                    << "In this case, this means the ansatz spaces of the cells as "
                    << "well as their derivatives. Each cell stores an index referring "
                    << "to their own ansatz spaces.",
        descriptionElement.setValue(description.view());
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const MeshData> rInstance) const {
        // Serialize basis functions.
        {
            GraphML::XMLElement element = rElement.addChild("basisFunctions");
            using SubSerializer = io::GraphML::Serializer<std::span<const Basis,polynomialOrder+1>>;
            SubSerializer subSerializer;
            subSerializer(element, rInstance.ansatzSpace().ansatzSet());
        }
    }
}; // struct GraphML::Serializer<MeshData>


/// @brief Deserializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Deserializer<MeshData>
    : public io::GraphML::DeserializerBase<MeshData>
{
    using io::GraphML::DeserializerBase<MeshData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair>) {
        // TODO
        (void)(pThis);
        CIE_THROW(NotImplementedException, elementName)
//        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
//
//        // Defer parsing to the array of ansatz spaces.
//        using SubDeserializer = io::GraphML::Deserializer<DynamicArray<Ansatz>>;
//        rThis.sax().push({
//            SubDeserializer::make(rThis._buffer, rThis.sax(), elementName),
//            SubDeserializer::onElementBegin,
//            SubDeserializer::onText,
//            SubDeserializer::onElementEnd
//        });
    }

    /// @brief This function is called when text block is parsed in the XML document.
    static void onText(void*,
                       std::string_view) {
        // No text data is expected for this class.
        CIE_THROW(Exception, "Unexpected text block while parsing mesh data.")
    }

    /// @brief This function is called when an element closing tag is parsed in the XML document.
    static void onElementEnd(void* pThis,
                             std::string_view elementName) {
        // TODO
        (void)(pThis);
        CIE_THROW(NotImplementedException, elementName)
//        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
//
//        // Move the parsed ansatz spaces from the buffer to the mesh data instance.
//        rThis.instance()._ansatzSpaces = std::move(rThis._buffer);
//
//        // Build derivatives from the parsed ansatz spaces.
//        std::transform(rThis.instance()._ansatzSpaces.begin(),
//                       rThis.instance()._ansatzSpaces.end(),
//                       std::back_inserter(rThis.instance()._ansatzDerivatives),
//                       [] (Ref<const Ansatz> rAnsatz) -> AnsatzDerivative {
//                            return rAnsatz.makeDerivative();
//                       });
//
//        // The parser's job is done => destroy it.
//        rThis.template release<Deserializer>(&rThis, elementName);
    }

private:
    DynamicArray<Ansatz> _buffer;
}; // struct GraphML::Deserializer<MeshData>


} // namespace cie::fem
