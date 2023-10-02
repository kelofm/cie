// Create a main() if CIE_DEFINE_TESTRUNNER_MAIN is defined
#ifdef CIE_DEFINE_TESTRUNNER_MAIN
#define CATCH_CONFIG_RUNNER
#endif

// Extend max error line size
//#define CATCH_CONFIG_CONSOLE_WIDTH 100

// --- External Includes ---
#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp> // Remove after catch2 v3 becomes available
#else
#include "catch2/catch_test_macros.hpp" // Uncomment after catch2 v3 becomes available
#include "catch2/catch_approx.hpp" // Uncomment after catch2 v3 becomes available
#include "catch2/catch_session.hpp" // Uncomment after catch2 v3 becomes available
namespace cie {
using Catch::Approx;
} // namespace cie
#endif

// --- Internal Includes ---
#include "packages/macros/inc/testing.hpp"
#include "packages/macros/inc/detail.hpp"



#define CIE_TEST_CASE(testCaseName, testCaseTag) TEST_CASE(testCaseName, testCaseTag)

#define CIE_TEST_CASE_METHOD(testCaseClass, testCaseName, testCaseTag) TEST_CASE_METHOD(testCaseClass, testCaseName, testCaseTag)

#define CIE_TEST_SECTION(sectionName) SECTION(sectionName)


// Test singles

#define CIE_TEST_CHECK(...) CIE_EXPAND_ARGUMENTS(CHECK(__VA_ARGS__))

#define CIE_TEST_CHECK_NOTHROW(...) CIE_EXPAND_ARGUMENTS(CHECK_NOTHROW(__VA_ARGS__))

#define CIE_TEST_CHECK_THROWS(...) CIE_EXPAND_ARGUMENTS(CHECK_THROWS(__VA_ARGS__))

#define CIE_TEST_CHECK_THROWS_AS(...) CIE_EXPAND_ARGUMENTS(CHECK_THROWS_AS(__VA_ARGS__))

#define CIE_TEST_REQUIRE(...) CIE_EXPAND_ARGUMENTS(REQUIRE(__VA_ARGS__))

#define CIE_TEST_REQUIRE_NOTHROW(...) CIE_EXPAND_ARGUMENTS(REQUIRE_NOTHROW(__VA_ARGS__))

#define CIE_TEST_REQUIRE_THROWS(...) CIE_EXPAND_ARGUMENTS(REQUIRE_THROWS(__VA_ARGS__))

#define CIE_TEST_REQUIRE_THROWS_AS(...) CIE_EXPAND_ARGUMENTS(REQUIRE_THROWS_AS(__VA_ARGS__))


#ifdef CIE_DEFINE_TESTRUNNER_MAIN
#ifdef CIE_TESTRUNNER_NAME

int main(int argc, char* argv[])
{
    auto p_logger = cie::utils::LoggerSingleton::getPtr(cie::getOutputPath()/ (std::string(CIE_STRINGIFY(CIE_TESTRUNNER_NAME)) + ".log"));
    const int result = Catch::Session().run(argc, argv);
    p_logger.reset();

    return result;
}
#endif
#endif