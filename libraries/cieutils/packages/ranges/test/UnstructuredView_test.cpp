// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/ranges/inc/UnstructuredView.hpp"


namespace cie::utils {


CIE_TEST_CASE("UnstructuredView", "[ranges]")
{
    CIE_TEST_CASE_INIT("UnstructuredView")

    int value = -2;
    const int constValue = -1;
    DynamicArray<int> array {0, 10};
    DynamicArray<int*> pointerArray {&array[0], &array[1]};
    DynamicArray<const int*> constPtrArray {&array[0], &array[1]};

    using TestView = UnstructuredView<int>;
    using ConstTestView = UnstructuredView<int,DynamicArray<const int*>>;

    CIE_TEST_CHECK_NOTHROW(UnstructuredView<int> {DynamicArray<int*> {&array[0], &value, &array[1]}});
    CIE_TEST_CHECK_NOTHROW(UnstructuredView<int,DynamicArray<const int*>> {DynamicArray<const int*> {&array[0], &value, &constValue, &array[1]}});

    {
        ConstTestView view;
        CIE_TEST_CHECK_NOTHROW(view = ConstTestView {pointerArray.begin(), pointerArray.end()});
        CIE_TEST_REQUIRE(view.size() == pointerArray.size());

        auto it = view.begin();
        for (int i=0; i<int(array.size()); ++i, ++it)
        {
            CIE_TEST_CHECK(view[i] == array[i]);
            CIE_TEST_CHECK(*(view.begin() + i) == array[i]);
        }
    }

    {
        ConstTestView view;
        CIE_TEST_CHECK_NOTHROW(view = ConstTestView {constPtrArray.begin(), constPtrArray.end()});
        CIE_TEST_REQUIRE(view.size() == pointerArray.size());

        auto it = view.begin();
        for (int i=0; i<int(array.size()); ++i, ++it)
        {
            CIE_TEST_CHECK(view[i] == array[i]);
            CIE_TEST_CHECK(*(view.begin() + i) == array[i]);
        }
    }

    {
        ConstTestView view;
        CIE_TEST_CHECK_NOTHROW(view = ConstTestView {array.begin(), array.end()});
        CIE_TEST_REQUIRE(view.size() == pointerArray.size());

        auto it = view.begin();
        for (int i=0; i<int(array.size()); ++i, ++it)
        {
            CIE_TEST_CHECK(view[i] == array[i]);
            CIE_TEST_CHECK(*(view.begin() + i) == array[i]);
        }
    }

    {
        TestView view;
        CIE_TEST_CHECK_NOTHROW(view = TestView {pointerArray.begin(), pointerArray.end()});
        CIE_TEST_REQUIRE(view.size() == pointerArray.size());

        auto it = view.begin();
        for (int i=0; i<int(array.size()); ++i, ++it)
        {
            CIE_TEST_CHECK(view[i] == array[i]);
            CIE_TEST_CHECK(*(view.begin() + i) == array[i]);

            view[i] = 10 - view[i];
        }

        CIE_TEST_CHECK(view[0] == 10);
        CIE_TEST_CHECK(view[1] == 0);
    }

    {
        TestView view;
        CIE_TEST_CHECK_NOTHROW(view = TestView {array.begin(), array.end()});
        CIE_TEST_REQUIRE(view.size() == pointerArray.size());

        auto it = view.begin();
        for (int i=0; i<int(array.size()); ++i, ++it)
        {
            CIE_TEST_CHECK(view[i] == array[i]);
            CIE_TEST_CHECK(*(view.begin() + i) == array[i]);

            view[i] = 10 * i;
        }

        CIE_TEST_CHECK(view[0] == 0);
        CIE_TEST_CHECK(view[1] == 10);
    }
} // CIE_TEST_CASE(UnstructuredView)


} // namespace cie::utils
