// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/testing/inc/ConstructorTracker.hpp"

// --- Internal Includes ---
#include "packages/stl_extension/inc/zip.hpp"

// --- STL Includes ---
#include <vector>
#include <array>
#include <map>
#include <list>


namespace cie::utils {


namespace testing {
template <class TIterator, class TDistance>
TIterator advanceAndReturn(TIterator it, TDistance distance)
{
    std::advance(it, distance);
    return it;
}
} // namespace testing


CIE_TEST_CASE("zip", "[stl_extension][!mayfail]")
{
    CIE_TEST_CASE_INIT("zip")

    {
        CIE_TEST_CASE_INIT("Containers of built-ins")

        {
            CIE_TEST_CASE_INIT("Uniform qualifiers")

            Size index = 0;
            DynamicArray<char> dynamicArray {'a', 'b', 'c'};
            std::vector<int> vector {0, 1, 2};
            std::array<float,3> array {0.0, 0.1, 0.2};
            std::map<char,float> map {{'a', 0}, {'b', 1}, {'c', 2}};
            std::list<double> list {0.0, 0.1, 0.2};

            for (const auto& [dynamicArrayItem, vectorItem, arrayItem, mapItem, listItem] : zip(dynamicArray, vector, array, map, list)) {
                CIE_TEST_CHECK(dynamicArrayItem == dynamicArray[index]);
                CIE_TEST_CHECK(vectorItem == vector[index]);
                CIE_TEST_CHECK(arrayItem == array[index]);
                CIE_TEST_CHECK(mapItem.first == testing::advanceAndReturn(map.begin(),index)->first);
                CIE_TEST_CHECK(mapItem.second == testing::advanceAndReturn(map.begin(),index)->second);
                CIE_TEST_CHECK(listItem == *testing::advanceAndReturn(list.begin(), index));
                ++index;
            }
        } // Uniform qualifiers

        {
            CIE_TEST_CASE_INIT("Mixed qualifiers")

            Size index = 0;
            DynamicArray<char> dynamicArray {'a', 'b', 'c'};
            const auto& constDynamicArray = dynamicArray;

            for (const auto& [c, cc] : zip(dynamicArray, constDynamicArray))
            {
                CIE_TEST_CHECK(c == dynamicArray[index]);
                CIE_TEST_CHECK(cc == dynamicArray[index]);
                ++index;
            }
        } // Mixed qualifiers

        {
            CIE_TEST_CASE_INIT("Mutable values")

            int index = 0;

            DynamicArray<int> dynamicArray {-1, -1, -1};
            std::vector<int> vector {-1, -1, -1};
            for (const auto& [arrayItem, vectorItem] : zip(dynamicArray, vector))
            {
                dynamicArray[index] = index;
                vector[index] = 10 - index;
                CIE_TEST_CHECK(arrayItem == index);
                CIE_TEST_CHECK(vectorItem == 10 - index);
                ++index;
            }

            // zip returns an rvalue, so its elements can only be referenced
            // with rvalue references
            for (auto&& [arrayItem, vectorItem] : zip(dynamicArray, vector))
            {
                arrayItem = 0;
                vectorItem = 0;
            }

            for (Size i=0; i<3; ++i)
            {
                CIE_TEST_CHECK(dynamicArray[i] == 0);
                CIE_TEST_CHECK(vector[i] == 0);
            }
        } // Mutable values

        {
            CIE_TEST_CASE_INIT("Mutable containers")

            int index = 0;
            DynamicArray<int> array0 {1, 2, 3};
            DynamicArray<int> array1 {4, 5, 6};

            // Make sure the iterators won't get invalidated
            array0.reserve(10);
            array1.reserve(10);

            auto zipped = zip(array0, array1);
            for (auto it=zipped.begin(); it!=zipped.end(); ++it) {
                auto tuple = *it;
                if (index < 7) {
                    // BTW don't ever do this
                    array0.push_back(index + 4);
                    array1.push_back(index + 7);
                }
                CIE_TEST_CHECK(std::get<0>(tuple) == index + 1);
                CIE_TEST_CHECK(std::get<1>(tuple) == index + 4);
                ++index;
            }
            CIE_TEST_CHECK(index == 10);
        } // Mutable containers

        {
            CIE_TEST_CASE_INIT("Mismatched sizes")

            int index = 0;
            DynamicArray<int> array3 {0, 1, 2};
            const DynamicArray<int> array10 {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

            for (const auto& [value0, value1] : zip(array3, array10))
            {
                CIE_TEST_CHECK(value0 == index);
                CIE_TEST_CHECK(value1 == index + 1);
                ++index;
            }
            CIE_TEST_CHECK(index == 3);

            index = 0;
            for (const auto& [value0, value1] : zip(array10, array3))
            {
                CIE_TEST_CHECK(value0 == index + 1);
                CIE_TEST_CHECK(value1 == index);
                ++index;
            }
            CIE_TEST_CHECK(index == 3);
        } // Mismatched sizes
    } // Containers of builtins

    {
        CIE_TEST_CASE_INIT("Move semantics")

        using Tracker = ConstructorTracker<>;
        DynamicArray<Tracker> array(3);
        for (auto&& [r_item, r_dummy] : zip(array, array))
        {
            CIE_TEST_CHECK(r_item.numberOfDefaultConstructorCalls() == 3);
            CIE_TEST_CHECK(r_item.numberOfActiveInstances() == 3);
            CIE_TEST_CHECK(r_item.numberOfMoveConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfCopyAssignmentCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfMoveConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfCopyConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfDestructorCalls() == 0);
            CIE_TEST_CHECK(r_dummy.numberOfDestructorCalls() == 0);
        }

        CIE_TEST_CHECK(Tracker::numberOfDefaultConstructorCalls() == 3);
        CIE_TEST_CHECK(Tracker::numberOfActiveInstances() == 3);
        CIE_TEST_CHECK(Tracker::numberOfMoveConstructorCalls() == 0);
        CIE_TEST_CHECK(Tracker::numberOfCopyAssignmentCalls() == 0);
        CIE_TEST_CHECK(Tracker::numberOfMoveConstructorCalls() == 0);
        CIE_TEST_CHECK(Tracker::numberOfCopyConstructorCalls() == 0);
        CIE_TEST_CHECK(Tracker::numberOfDestructorCalls() == 0);

        array.clear();
        Tracker::reset();

        /// @todo prevent copying on moved or rvalue containers
        for (auto&& [r_item] : zip(DynamicArray<Tracker>(3))) {
            CIE_TEST_CHECK(r_item.numberOfDefaultConstructorCalls() == 3);
            CIE_TEST_CHECK(r_item.numberOfActiveInstances() == 3);
            CIE_TEST_CHECK(r_item.numberOfMoveConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfCopyAssignmentCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfMoveConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfCopyConstructorCalls() == 0);
            CIE_TEST_CHECK(r_item.numberOfDestructorCalls() == 0);
        }
    }
}


} // namespace cie::utils