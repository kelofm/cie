#pragma once

// --- Utility Includes ---
#include "packages/types/inc/tags.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"

// --- STL Includes ---
#include <vector>
#include <span>


namespace cie::linalg {


template <class T, TagLike TTag = tags::SMP>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
class DefaultSpace {
public:
    using Value = T;

    using Vector = std::vector<T>;

    using VectorView = std::span<T>;

    using ConstVectorView = std::span<const T>;

    constexpr DefaultSpace() noexcept = default;

    DefaultSpace(Ref<mp::ThreadPoolBase> rThreads)
    requires std::is_same_v<TTag,tags::SMP>;

    static VectorView view(Ref<Vector> rVector) noexcept;

    static ConstVectorView view(Ref<const Vector> rVector) noexcept;

    static std::size_t size(ConstVectorView view) noexcept;

    static Vector makeVector(std::size_t size);

    Value innerProduct(ConstVectorView left, ConstVectorView right) const;

    void scale(VectorView target, ConstVectorView source, Value scale) const;

    void scale(VectorView view, Value value) const;

    void add(VectorView target, ConstVectorView source, Value scale) const;

    void assign(VectorView target, ConstVectorView source) const;

    void fill(VectorView view, Value value) const;

private:
    OptionalRef<mp::ThreadPoolBase> _maybeThreads;
}; // class DefaultSpace


} // namespace cie::linalg
