#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <span>


namespace cie {


enum class ReorderingStrategy {
    None,
    CuthillMcKee,
    ReverseCuthillMcKee
}; // enum class ReorderingStrategy


template <concepts::Integer TIndex, concepts::Numeric TValue>
void makeReordering(
    std::span<TIndex> map,
    std::span<const TIndex> rowExtents,
    std::span<const TIndex> columnIndices,
    std::span<const TValue> entries,
    ReorderingStrategy strategy = ReorderingStrategy::CuthillMcKee,
    OptionalRef<mp::ThreadPoolBase> rMaybePool = {});


template <concepts::Integer TIndex>
void reverseReorder(
    std::span<TIndex> map,
    std::span<TIndex> buffer,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads = {});


template <concepts::Integer TIndex, concepts::Numeric TValue, bool Reverse = false>
void reorder(
    std::span<const TIndex> map,
    std::span<TIndex> rowExtents,
    std::span<TIndex> columnIndices,
    std::span<TValue> entries,
    OptionalRef<mp::ThreadPoolBase> rMaybePool = {});


template <concepts::Integer TIndex, concepts::Numeric TValue>
void reverseReorder(
    std::span<const TIndex> map,
    std::span<TIndex> rowExtents,
    std::span<TIndex> columnIndices,
    std::span<TValue> entries,
    OptionalRef<mp::ThreadPoolBase> rMaybePool = {});


template <concepts::Integer TIndex, concepts::Numeric TValue, bool Reverse = false>
void reorder(
    std::span<const TIndex> map,
    std::span<TValue> array,
    OptionalRef<mp::ThreadPoolBase> rMaybePool = {});


template <concepts::Integer TIndex, concepts::Numeric TValue>
void reverseReorder(
    std::span<const TIndex> map,
    std::span<TValue> array,
    OptionalRef<mp::ThreadPoolBase> rMaybePool = {});


} // namespace cie
