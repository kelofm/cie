// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/concepts/inc/partial_specialization.hpp"


namespace cie::concepts {


CIE_DEFINE_INVALID_CLASS_TEMPLATE_TO_SPECIALIZE( PartialSpecializationTestClass )

template <Integer T>
class PartialSpecializationTestClass<T>
{
public:
    void operator()() {}
};


CIE_TEST_CASE( "Partial Specialization", "[concepts]" )
{
    CIE_TEST_CASE_INIT( "Partial Specialization" )

    // This should compile
    PartialSpecializationTestClass<int> test;
    test();

    // This should not compile
    //PartialSpecializationTestClass<double> fail;
}


} // namespace cie::utils