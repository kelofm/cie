// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/stl_extension/inc/StackArray.hpp"


namespace cie {


CIE_TEST_CASE("StackArray", "[stl_extension]")
{
    CIE_TEST_CASE_INIT("StackArray")

    CIE_TEST_SECTION("std::size_t") {
        CIE_TEST_CASE_INIT("std::size_t")
        constexpr std::size_t maxSize = 10;
        using T = StackArray<std::size_t,maxSize>;

        CIE_TEST_REQUIRE_NOTHROW(T());
        T array;
        CIE_TEST_CHECK(array.capacity() == 10ul);
        CIE_TEST_REQUIRE(array.size() == 0ul);
        CIE_TEST_CHECK(array.empty());
        CIE_TEST_CHECK(array.begin() == array.end());
        CIE_TEST_CHECK(array.cbegin() == array.cend());
        CIE_TEST_CHECK(array.rbegin() == array.rend());
        CIE_TEST_CHECK(array.crbegin() == array.crend());

        for (std::size_t count=0ul; count<=maxSize; ++count) {
            bool result = false;

            CIE_TEST_CHECK_NOTHROW(result = array.reserve(count));
            CIE_TEST_REQUIRE(result);

            CIE_TEST_CHECK_NOTHROW(result = array.resize(count));
            CIE_TEST_REQUIRE(result);
            CIE_TEST_CHECK(array.size() == count);
            CIE_TEST_CHECK(std::distance(array.begin(), array.end()) == static_cast<std::ptrdiff_t>(count));

            for (std::size_t iItem=0ul; iItem<count; ++iItem) {
                CIE_TEST_CHECK_NOTHROW(array.at(iItem) = iItem * iItem);
            }

            for (std::size_t iItem=0ul; iItem<count; ++iItem) {
                CIE_TEST_CHECK(array[iItem] == iItem * iItem);
            }
        }

        CIE_TEST_CHECK(!array.reserve(maxSize + 1));
        CIE_TEST_CHECK(!array.reserve(maxSize + 1e3));

        CIE_TEST_CHECK_NOTHROW(array.clear());
        CIE_TEST_CHECK(array.capacity() == 10ul);
        CIE_TEST_REQUIRE(array.size() == 0ul);
        CIE_TEST_CHECK(array.empty());
        CIE_TEST_CHECK(array.begin() == array.end());
        CIE_TEST_CHECK(array.cbegin() == array.cend());
        CIE_TEST_CHECK(array.rbegin() == array.rend());
        CIE_TEST_CHECK(array.crbegin() == array.crend());

        for (std::size_t iItem=0ul; iItem<maxSize; ++iItem) {
            if (iItem % 2) {
                CIE_TEST_CHECK(array.push_back(iItem * iItem));
            } else {
                CIE_TEST_CHECK(array.emplace_back(iItem * iItem));
            }
            CIE_TEST_CHECK(array.back() == iItem * iItem);
        }

        for (std::size_t iItem=0ul; iItem<maxSize; ++iItem) {
            CIE_TEST_CHECK_NOTHROW(array.at(iItem) = iItem * iItem);
        }
    }
}


} // namespace cie
