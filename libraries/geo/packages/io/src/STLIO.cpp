// --- GEO Includes ---
#include "packages/io/inc/STLIO.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <array>
#include <cstdint>
#include <algorithm>
#include <format>


namespace cie::io {


template <class TV, unsigned D>
struct STLIO::Input<TV,D>::Impl {
    Ptr<std::istream> pStream;

    std::size_t triangleCount;
};


template <class TV, unsigned D>
STLIO::Input<TV,D>::Input(Ref<std::istream> rStream)
    : _pImpl(new Impl {
        .pStream        = &rStream,
        .triangleCount  = 0ul}) {
            CIE_BEGIN_EXCEPTION_TRACING
                // Discard header (80 bytes).
                rStream.ignore(80);

                // Read number of triangles.
                std::uint32_t triangleCount = 0;
                rStream.read(
                    reinterpret_cast<char*>(&triangleCount),
                    sizeof(triangleCount));
                _pImpl->triangleCount = triangleCount;
            CIE_END_EXCEPTION_TRACING
}


template <class TV, unsigned D>
STLIO::Input<TV,D>::~Input() = default;


template <class TV, unsigned D>
std::size_t STLIO::Input<TV,D>::triangleCount() const noexcept {
    return _pImpl->triangleCount;
}


template <class TV, unsigned D, bool ReadNormals>
void readSTL(
    Ref<std::istream> rStream,
    std::span<TV> coordinates,
    [[maybe_unused]] std::span<TV> normals,
    std::size_t triangleCount) {
        Ptr<TV> itCoordinate = coordinates.data();
        [[maybe_unused]] Ptr<TV> itNormal = normals.data();
        std::array<float,3> buffer;
        for (std::size_t iTriangle=0ul; iTriangle<triangleCount; ++iTriangle) {
            if constexpr (ReadNormals) {
                rStream.read(
                    reinterpret_cast<Ptr<char>>(buffer.data()),
                    3 * sizeof(float));
                std::copy_n(
                    buffer.data(),
                    3,
                    itNormal);
                itNormal += 3;
            } else {
                rStream.ignore(3 * sizeof(float));
            }

            for (std::size_t iVertex=0ul; iVertex<3; ++iVertex) {
                rStream.read(
                    reinterpret_cast<Ptr<char>>(buffer.data()),
                    D * sizeof(float));
                std::copy_n(
                    buffer.data(),
                    D,
                    itCoordinate);
                if (D != 3) rStream.ignore((3 - D) * sizeof(float));
                itCoordinate += D;
            }

            rStream.ignore(sizeof(std::uint16_t));
        }
}


template <class TV, unsigned D>
void STLIO::Input<TV,D>::execute(
    std::span<TV> coordinates,
    std::span<TV> normals) {
        CIE_CHECK(
            coordinates.size() == 3 * D * _pImpl->triangleCount,
            std::format(
                "expecting coordinate array of size {} for {} {}D triangles, but got {}",
                3 * D * _pImpl->triangleCount,
                _pImpl->triangleCount,
                D,
                coordinates.size()))
        CIE_CHECK(
            normals.empty() || normals.size() == 3 * _pImpl->triangleCount,
            std::format(
                "expecting normal array of size {} for {} {}D triangles, but got {}",
                3 * _pImpl->triangleCount,
                _pImpl->triangleCount,
                D,
                normals.size()))

        CIE_BEGIN_EXCEPTION_TRACING
            if (normals.empty()) {
                readSTL<TV,D,false>(
                    *_pImpl->pStream,
                    coordinates,
                    normals,
                    _pImpl->triangleCount);
            } else {
                readSTL<TV,D,true>(
                    *_pImpl->pStream,
                    coordinates,
                    normals,
                    _pImpl->triangleCount);
            }
        CIE_END_EXCEPTION_TRACING
}


template class STLIO::Input<float,2>;
template class STLIO::Input<float,3>;
template class STLIO::Input<double,2>;
template class STLIO::Input<double,3>;


template <class TV, unsigned D>
struct STLIO::Output<TV,D>::Impl {
    Ptr<std::ostream> pStream;
}; // struct STLIO::Output::Impl


template <class TV, unsigned D>
STLIO::Output<TV,D>::Output(Ref<std::ostream> rStream)
    : _pImpl(new Impl {.pStream = &rStream})
{}


template <class TV, unsigned D>
STLIO::Output<TV,D>::~Output() = default;


template <class TV, unsigned D, bool WriteNormals>
void writeSTL(
    Ref<std::ostream> rStream,
    std::span<const TV> coordinates,
    [[maybe_unused]] std::span<const TV> normals) {
        // Write header.
        {
            std::array<char,80> header;
            std::fill(header.begin(), header.end(), '0');
            rStream.write(header.data(), header.size());
        }

        // Write number of triangles.
        const std::uint32_t triangleCount = coordinates.size() / D / 3;
        rStream.write(
            reinterpret_cast<const char*>(&triangleCount),
            sizeof(std::uint32_t));

        // Write triangles.
        for (std::uint32_t iTriangle=0ul; iTriangle<triangleCount; ++iTriangle) {
            // Collect and write the normal.
            std::array<float,3> normal;
            if constexpr (WriteNormals) {
                std::copy_n(
                    normals.begin() + iTriangle * 3,
                    3,
                    normal.begin());
            } else {
                normal[0] = 0.0f;
                normal[1] = 0.0f;
                normal[2] = 1.0f;
            }
            rStream.write(
                reinterpret_cast<Ptr<const char>>(normal.data()),
                sizeof(float) * normal.size());

            // Collect and write vertices.
            std::array<float,9> vertices;
            const std::size_t iCoordinateBegin = iTriangle * 3 * D;
            if constexpr (D == 3) {
                std::copy_n(
                    coordinates.begin() + iCoordinateBegin,
                    9,
                    vertices.begin());
            } else {
                for (unsigned iVertex=0u; iVertex<3; ++iVertex) {
                    for (unsigned iDimension=0u; iDimension<D; ++iDimension)
                        vertices[3 * iVertex + iDimension] = static_cast<float>(coordinates[iCoordinateBegin + D * iVertex + iDimension]);
                    for (unsigned iDimension=D; iDimension<3; ++iDimension)
                        vertices[3 * iVertex + iDimension] = 0.0f;
                }
            }
            rStream.write(
                reinterpret_cast<Ptr<const char>>(vertices.data()),
                sizeof(float) * vertices.size());

            // Write attribute byte count.
            const std::uint16_t attributeByteCount = 1;
            rStream.write(
                reinterpret_cast<Ptr<const char>>(&attributeByteCount),
                sizeof(attributeByteCount));
        } // for iTriangle in range(triangleCount)
}


template <class TV, unsigned D>
void STLIO::Output<TV,D>::execute(
    std::span<const TV> coordinates,
    std::span<const TV> normals) {
        CIE_BEGIN_EXCEPTION_TRACING
            if (normals.empty())
                writeSTL<TV,D,false>(*_pImpl->pStream, coordinates, normals);
            else
                writeSTL<TV,D,true>(*_pImpl->pStream, coordinates, normals);
        CIE_END_EXCEPTION_TRACING
}


template class STLIO::Output<float,2>;
template class STLIO::Output<float,3>;
template class STLIO::Output<double,2>;
template class STLIO::Output<double,3>;


} // namespace cie::io
