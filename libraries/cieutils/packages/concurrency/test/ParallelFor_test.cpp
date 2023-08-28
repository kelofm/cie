// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"

// --- STL Includes ---
#include <vector>
#include <string>
#include <numeric>


namespace cie::mp {


CIE_TEST_CASE("ParallelFor", "[concurrency]")
{
    CIE_TEST_CASE_INIT("ParallelFor")

    const Size indexMin       = 5;
    const Size indexMax       = 1005;
    const Size stepSize       = 2;

    {
        CIE_TEST_CASE_INIT("index loop - empty storage")

        DynamicArray<int> array(2 * indexMax, 0);

        ParallelFor<>()(
            indexMin,
            indexMax,
            stepSize,
            [&array](Size index)
            {++array[index];}
        );

        for (Size i=0; i<array.size(); ++i)
        {
            const bool touched = indexMin <= i && i < indexMax && (i - indexMin) % stepSize == 0;
            if (touched)
                CIE_TEST_CHECK(array[i] == 1);
            else
                CIE_TEST_CHECK(array[i] == 0);
        }
    }

    {
        CIE_TEST_CASE_INIT("index loop")

        DynamicArray<int> array(2 * indexMax, 1);
        int firstPrivateVar = 0;

        auto results = ParallelFor<>::firstPrivate(firstPrivateVar)(
            indexMin,
            indexMax,
            stepSize,
            [&array](Size index, int& privateVar)
            {
                privateVar += array[index];
                array[index] = -1;
            }
        ).getPool().getStorage();

        const int sum = std::accumulate(
            results.begin(),
            results.end(),
            0,
            [](int lhs, const auto& r_storage) {return lhs + r_storage.template get<0>();}
        );

        CIE_TEST_CHECK(sum == (indexMax - indexMin) / stepSize);

        for (Size i=0; i<array.size(); ++i)
        {
            const bool touched = indexMin <= i && i < indexMax && (i - indexMin) % stepSize == 0;
            if (touched)
                CIE_TEST_CHECK(array[i] == -1);
            else
                CIE_TEST_CHECK(array[i] == 1);
        }
    }

    {
        CIE_TEST_CASE_INIT("object loop")

        std::vector<std::string> objects;
        objects.reserve(indexMax - indexMin);

        for (Size i=0; i<indexMax-indexMin; ++i)
            objects.push_back(std::to_string(i));

        std::string suffix = "_";

        ParallelFor<>::firstPrivate(suffix)(
            objects,
            [](std::string& r_item, std::string& r_suffix) -> void
            {
                r_item += r_suffix;
            }
        );

        for (Size i=0; i<objects.size(); ++i)
            CIE_TEST_CHECK(
                objects[i]
                ==
                std::to_string(i) + "_"
            );
    }
}


} // namespace cie::mp