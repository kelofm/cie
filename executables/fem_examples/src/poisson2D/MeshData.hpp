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

    template <class T>
    static sycl::usm_allocator<T,sycl::usm::alloc::device> makeDeviceAllocator() {
        return sycl::usm_allocator<T,sycl::usm::alloc::device>(SYCLSingleton::getQueue());
    }

private:
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
}; // class SYCLSingleton


std::optional<SYCLSingleton::Impl> SYCLSingleton::_maybeImpl = {};


/// @brief Data structure common to the entire @ref Graph "mesh".
class MeshData {
public:
    MeshData()
        : _polynomialCoefficients(SYCLSingleton::makeSharedAllocator<Scalar>()),
          _basisFunctions(SYCLSingleton::makeSharedAllocator<Basis>()),
          _basisDerivatives(SYCLSingleton::makeSharedAllocator<Basis::Derivative>()),
          _ansatzSpaces(SYCLSingleton::makeSharedAllocator<Ansatz>()),
          _ansatzDerivatives(SYCLSingleton::makeSharedAllocator<AnsatzDerivative>()),
          _quadraturePointSets(),
          _buffer()
    {}

    MeshData(std::span<unsigned> integrationOrders,
             RightRef<DynamicSharedArray<Scalar>> rPolynomialCoefficients,
             RightRef<DynamicSharedArray<Basis>> rBasisFunctions,
             RightRef<DynamicSharedArray<Basis::Derivative>> rBasisDerivatives,
             RightRef<DynamicSharedArray<Ansatz>> rAnsatzSpaces,
             RightRef<DynamicSharedArray<AnsatzDerivative>> rAnsatzDerivatives)
        : _polynomialCoefficients(std::move(rPolynomialCoefficients)),
          _basisFunctions(std::move(rBasisFunctions)),
          _basisDerivatives(std::move(rBasisDerivatives)),
          _ansatzSpaces(std::move(rAnsatzSpaces)),
          _ansatzDerivatives(std::move(rAnsatzDerivatives)),
          _quadraturePointSets()
    {
        // Cache quadrature points.
        _quadraturePointSets.resize(integrationOrders.size());
        CellData dummyElement;
        for (std::size_t iSet=0u; iSet<integrationOrders.size(); ++iSet) {
            const std::size_t integrationOrder = integrationOrders[iSet];

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
            auto& rSet = _quadraturePointSets[iSet];
            rSet.resize(initialSize);

            while (true) {
                // Extend the container if necessary.
                if (!(pointCount < rSet.size())) {
                    rSet.resize(growFactor * rSet.size());
                }

                std::span<QuadraturePoint<Dimension,Scalar>> targetSpan(
                    rSet.data() + pointCount,
                    rSet.data() + rSet.size());

                // Request a new batch of quadrature points.
                const auto newPointCount = generator.generate(dummyElement, targetSpan);
                pointCount += newPointCount;

                if (!newPointCount) break;
            }

            rSet.resize(pointCount);
        } // for iSet in range(integrationOrders.size())

        this->fitBuffer();
    }

    std::span<const Scalar> polynomialCoefficients() const noexcept {
        return _polynomialCoefficients;
    }

    std::span<Scalar> polynomialCoefficients() noexcept {
        return _polynomialCoefficients;
    }

    std::span<const Ansatz> ansatzSpaces() const noexcept {
        return _ansatzSpaces;
    }

    std::span<const AnsatzDerivative> ansatzDerivatives() const noexcept {
        return _ansatzDerivatives;
    }

    CachedQuadraturePointFactory<Dimension,Scalar> makeQuadratureRule(std::size_t iAnsatz) const {
        return CachedQuadraturePointFactory<Dimension,Scalar>(
            std::span<const QuadraturePoint<Dimension,Scalar>>(
                _quadraturePointSets[iAnsatz].data(),
                _quadraturePointSets[iAnsatz].size())
        );
    }

private:
    void fitBuffer() {
        std::size_t bufferSize = 0ul;
        for (Ref<const Ansatz> rAnsatz : _ansatzSpaces)
            bufferSize += rAnsatz.getMinBufferSize();
        for (Ref<const AnsatzDerivative> rDerivative : _ansatzDerivatives)
            bufferSize += rDerivative.getMinBufferSize();

        _buffer.resize(bufferSize);

        bufferSize = 0ul;
        for (Ref<Ansatz> rAnsatz : _ansatzSpaces) {
            const std::size_t partialSize = rAnsatz.getMinBufferSize();
            rAnsatz.setBuffer({
                _buffer.data() + bufferSize,
                partialSize});
            bufferSize += partialSize;
        }
        for (Ref<AnsatzDerivative> rDerivative : _ansatzDerivatives) {
            const std::size_t partialSize = rDerivative.getMinBufferSize();
            rDerivative.setBuffer({
                _buffer.data() + bufferSize,
                partialSize});
            bufferSize += partialSize;
        }
    }

    friend struct io::GraphML::Serializer<MeshData>;

    friend struct io::GraphML::Deserializer<MeshData>;

    DynamicSharedArray<Scalar> _polynomialCoefficients;

    DynamicSharedArray<Basis> _basisFunctions;

    DynamicSharedArray<Basis::Derivative> _basisDerivatives;

    /// @brief Collection of all ansatz spaces the contained cells can refer to.
    DynamicSharedArray<Ansatz> _ansatzSpaces;

    /// @brief Collection of all ansatz spaces' derivatives the contained cells can refer to.
    DynamicSharedArray<AnsatzDerivative> _ansatzDerivatives;

    /// @brief Sets of quadrature points for a default local hypercube.
    /// @details These quadrature points are used while constructing
    ///          cell-specific quadrature rules.
    DynamicArray<DynamicArray<QuadraturePoint<Dimension,Scalar>>> _quadraturePointSets;

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
        // Serialize ansatz spaces.
        {
            GraphML::XMLElement containerElement = rElement.addChild("ansatzSpaces");
            for (Ref<const Ansatz> rAnsatz : rInstance._ansatzSpaces) {
                GraphML::XMLElement ansatzElement = containerElement.addChild("ansatzSpace");
                using SubSerializer = io::GraphML::Serializer<std::span<const Basis>>;
                SubSerializer subSerializer;
                subSerializer(ansatzElement, rAnsatz.ansatzSet());
            } // for rAnsatz in rInstance._ansatzSpaces
        }

        // Serialize ansatz derivatives.
        {
            GraphML::XMLElement containerElement = rElement.addChild("ansatzDerivatives");
            for (Ref<const Ansatz> rAnsatz : rInstance._ansatzSpaces) {
                {
                    GraphML::XMLElement element = containerElement.addChild("basisFunctions");
                    using SubSerializer = io::GraphML::Serializer<std::span<const Basis>>;
                    SubSerializer subSerializer;
                    subSerializer(element, rAnsatz.ansatzSet());
                }
                {
                    GraphML::XMLElement element = containerElement.addChild("basisDerivatives");
                    using SubSerializer = io::GraphML::Serializer<std::span<const Basis>>;
                    SubSerializer subSerializer;
                    subSerializer(element, rAnsatz.ansatzSet());
                }
            } // for rAnsatz in rInstance._ansatzSpaces
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
