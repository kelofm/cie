// --- FEM Includes ---
#include "packages/utilities/inc/AttributeContainer.hpp"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"


namespace cie::fem {


CIE_TEST_CASE("AttributeContainer", "[container]")
{
    CIE_TEST_CASE_INIT("AttributeContainer")

    AttributeContainer<double,int> attributes;
    CIE_TEST_CHECK(attributes.empty());
    CIE_TEST_CHECK(attributes.size() == 0);
    attributes.resize(4);

    CIE_TEST_CHECK(attributes.size() == 4);
    CIE_TEST_CHECK(!attributes.empty());

    for (Size i=0; i<attributes.size(); ++i) {
        attributes.at<double>(i) = 1.0 + i * 0.5;
        attributes.at<int>(i) = i + 1;

        CIE_TEST_CHECK(attributes.at<double>(i) == Approx(1.0 + i * 0.5));
        CIE_TEST_CHECK(attributes.at<int>(i) == i + 1);
    }

    const auto visitor = [](double d, int i) {
        return d * i;
    };
    for (Size i=0; i<attributes.size(); ++i) {
        CIE_TEST_CHECK(attributes.visit<double,int>(i, visitor) == Approx((1.0 + i * 0.5) * (i + 1)));
    }

    CIE_TEST_CHECK_NOTHROW(attributes.erase(1, 3));
    CIE_TEST_CHECK(!attributes.empty());
    CIE_TEST_CHECK(attributes.size() == 2);
    CIE_TEST_CHECK(attributes.at<double>(0) == Approx(1.0 + 0 * 0.5));
    CIE_TEST_CHECK(attributes.at<int>(0) == 0 + 1);
    CIE_TEST_CHECK(attributes.at<double>(1) == Approx(1.0 + 3 * 0.5));
    CIE_TEST_CHECK(attributes.at<int>(1) == 3 + 1);

    CIE_TEST_CHECK_NOTHROW(attributes.erase(0));
    CIE_TEST_CHECK(!attributes.empty());
    CIE_TEST_CHECK(attributes.size() == 1);
    CIE_TEST_CHECK(attributes.at<double>(0) == Approx(1.0 + 3 * 0.5));
    CIE_TEST_CHECK(attributes.at<int>(0) == 3 + 1);

    CIE_TEST_CHECK_NOTHROW(attributes.clear());
    CIE_TEST_CHECK(attributes.empty());
    CIE_TEST_CHECK(attributes.size() == 0);

    const auto oldCapacity = attributes.capacity();
    attributes.reserve(oldCapacity * 2);
    CIE_TEST_CHECK(attributes.capacity() == oldCapacity * 2);
    CIE_TEST_CHECK(attributes.empty());
    CIE_TEST_CHECK(attributes.size() == 0);

    attributes.push_back(12.0, 144);
    CIE_TEST_CHECK(!attributes.empty());
    CIE_TEST_REQUIRE(attributes.size() == 1);
    CIE_TEST_CHECK(attributes.at<double>(0) == 12.0);
    CIE_TEST_CHECK(attributes.at<int>(0) == 144);

    const std::tuple<double,int> item = attributes.get(0);
    CIE_TEST_CHECK(std::get<0>(item) == Approx(12.0));
    CIE_TEST_CHECK(std::get<1>(item) == 144);
}


} // namespace cie::fem
