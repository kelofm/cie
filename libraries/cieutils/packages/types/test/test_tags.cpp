// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/types/inc/tags.hpp"


namespace cie::tags {


CIE_TEST_CASE("Tags", "[types]")
{
    CIE_TEST_CHECK(Null::isCompatibleWith<Null,Eager,Lazy,Serial,SMP,MPI>());

    CIE_TEST_CHECK(!Eager::isCompatibleWith<Lazy>());
    CIE_TEST_CHECK(Eager::isCompatibleWith<Null,Eager,Serial,SMP,MPI>());

    CIE_TEST_CHECK(!Lazy::isCompatibleWith<Eager>());
    CIE_TEST_CHECK(Lazy::isCompatibleWith<Null,Lazy,Serial,SMP,MPI>());

    CIE_TEST_CHECK((!Serial::isCompatibleWith<SMP>() && !Serial::isCompatibleWith<MPI>()));
    CIE_TEST_CHECK(Serial::isCompatibleWith<Null,Eager,Lazy,Serial>());

    CIE_TEST_CHECK((!SMP::isCompatibleWith<Serial>() && !SMP::isCompatibleWith<MPI>()));
    CIE_TEST_CHECK(SMP::isCompatibleWith<Null,Eager,Lazy,SMP>());

    CIE_TEST_CHECK((!MPI::isCompatibleWith<Serial>() && !MPI::isCompatibleWith<SMP>()));
    CIE_TEST_CHECK(MPI::isCompatibleWith<Null,Eager,Lazy,MPI>());
}


} // namespace cie::tags
