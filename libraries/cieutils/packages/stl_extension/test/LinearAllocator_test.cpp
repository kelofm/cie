// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/LinearAllocator.hpp"

// --- STL Includes ---
#include <vector>
#include <numeric>


namespace cie {


struct LinearAllocatorTestItem {
    std::size_t big;

    std::byte small;

    constexpr LinearAllocatorTestItem() noexcept
        : LinearAllocatorTestItem(0ul)
    {}

    constexpr LinearAllocatorTestItem(std::size_t i) noexcept
        :   big(i),
            small()
    {}
}; // struct LinearAllocatorTestItem


CIE_TEST_CASE("LinearAllocator", "[stl_extension]") {
    CIE_TEST_CASE_INIT("LinearAllocator")

    #define CIE_TEST_LINEAR_ALLOCATOR(T)                                        \
        std::vector<T> buffer(5);                                               \
        LinearArena arena {                                                     \
            .buffer = {                                                         \
                reinterpret_cast<std::byte*>(buffer.data()),                    \
                buffer.size() * sizeof(T)},                                     \
            .state = 0ul};                                                      \
        LinearAllocator<T> allocator(arena);                                    \
        std::vector<T,LinearAllocator<T>> linearBuffer(std::move(allocator));   \
        CIE_TEST_REQUIRE_NOTHROW(linearBuffer.resize(buffer.size()));           \
        std::iota(linearBuffer.begin(), linearBuffer.end(), 1);

    {
        CIE_TEST_CASE_INIT("unsigned short")
        CIE_TEST_LINEAR_ALLOCATOR(unsigned short)
        for (unsigned iItem=0u; iItem<static_cast<unsigned short>(buffer.size()); ++iItem)
            CIE_TEST_CHECK(buffer[iItem] == iItem + 1);
    }

    {
        CIE_TEST_CASE_INIT("int")
        CIE_TEST_LINEAR_ALLOCATOR(int)
        for (int iItem=0u; iItem<static_cast<int>(buffer.size()); ++iItem)
            CIE_TEST_CHECK(buffer[iItem] == iItem + 1);
    }

    {
        CIE_TEST_CASE_INIT("long long")
        CIE_TEST_LINEAR_ALLOCATOR(long long)
        for (long long iItem=0u; iItem<static_cast<long long>(buffer.size()); ++iItem)
            CIE_TEST_CHECK(buffer[iItem] == iItem + 1);
    }

    {
        CIE_TEST_CASE_INIT("float")
        CIE_TEST_LINEAR_ALLOCATOR(float)
        for (std::size_t iItem=0u; iItem<buffer.size(); ++iItem)
            CIE_TEST_CHECK(buffer[iItem] == iItem + 1);
    }

    #undef CIE_TEST_LINEAR_ALLOCATOR

} // "LinearAllocator"


} // namespace cie
