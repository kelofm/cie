// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Range.hpp"


namespace cie::ct {


template <class ...Ts>
struct RangeTest {};


CIE_TEST_CASE("Range", "[compile_time]")
{

    CIE_TEST_CASE_INIT("Range")

    using Full = RangeTest<void,char,char,int,double,RangeTest<>>;
    using RangeBase = Range<void,char,char,int,double,RangeTest<>>;

    using Ref_Begin_3 = RangeTest<void,char,char>;
    using Ref_3_End = RangeTest<int,double,RangeTest<>>;
    using Ref_1_3 = RangeTest<char,char>;
    using Ref_First = RangeTest<void>;
    using Ref_Last = RangeTest<RangeTest<>>;
    using Ref_Full = Full;
    using Ref_Empty = RangeTest<>;

    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<0,3>::Apply<RangeTest>, Ref_Begin_3>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<3,6>::Apply<RangeTest>, Ref_3_End>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<1,3>::Apply<RangeTest>, Ref_1_3>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<0,1>::Apply<RangeTest>, Ref_First>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<5,6>::Apply<RangeTest>, Ref_Last>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<0,6>::Apply<RangeTest>, Ref_Full>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<0,0>::Apply<RangeTest>, Ref_Empty>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<3,3>::Apply<RangeTest>, Ref_Empty>);
    CIE_TEST_CHECK(std::is_same_v<RangeBase::Subrange<5,5>::Apply<RangeTest>, Ref_Empty>);
}


} // namespace cie::ct