#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <cstddef>
#include <span>


namespace cie::io {


/// @brief Concept for streams compatible with binary deserialization.
/// @details A type satisfies this concept if it exposes a read function callable as
///          @code{.cpp} rStream.read(Ptr<char> pData, std::size_t byteCount) @endcode.
///          Both @p std::istream and @ref ByteDeserializer satisfy it.
/// @ingroup cieutils
template <class TStream>
concept BinaryDeserializerStreamLike
= requires (Ref<TStream> rStream, Ptr<char> pData, std::size_t byteCount) {
    {rStream.read(pData, byteCount)};
}; // concept BinaryDeserializerStreamLike


/// @brief Deserializer over a contiguous range of bytes.
/// @details Unlike @p std::istream , this class is trivially copyable and holds no
///          state beyond a triplet of pointers, which makes it usable in SYCL kernels.
///          Reading past the end of the range terminates the program.
/// @ingroup cieutils
class ByteDeserializer {
public:
    ByteDeserializer() noexcept = default;

    ByteDeserializer(
        Ptr<const std::byte> pData,
        std::size_t byteCount) noexcept
            :   _pBegin(pData),
                _pEnd(pData + byteCount),
                _pCurrent(pData)
    {}

    ByteDeserializer(std::span<const std::byte> rBytes) noexcept
        : ByteDeserializer(rBytes.data(), rBytes.size())
    {}

    void read(
        Ptr<char> pData,
        std::size_t byteCount) noexcept {
            if (byteCount == 0ul) return;
            if (static_cast<std::ptrdiff_t>(byteCount) <= _pEnd - _pCurrent) {
                for (std::size_t iByte=0ul; iByte<byteCount; ++iByte)
                    pData[iByte] = reinterpret_cast<Ptr<const char>>(_pCurrent)[iByte];
                _pCurrent += byteCount;
            } else {
                __builtin_trap();
            }
    }

    [[nodiscard]] std::size_t remaining() const noexcept {
        if (_pCurrent == nullptr) return 0ul;
        return static_cast<std::size_t>(_pEnd - _pCurrent);
    }

    [[nodiscard]] std::size_t consumed() const noexcept {
        if (_pCurrent == nullptr) return 0ul;
        return static_cast<std::size_t>(_pCurrent - _pBegin);
    }

private:
    Ptr<const std::byte> _pBegin = nullptr;

    Ptr<const std::byte> _pEnd = nullptr;

    Ptr<const std::byte> _pCurrent = nullptr;
}; // class ByteDeserializer


static_assert(BinaryDeserializerStreamLike<ByteDeserializer>);


} // namespace cie::io
