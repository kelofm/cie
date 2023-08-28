// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "cmake_variables.hpp"


namespace cie::utils {


CIE_TEST_CASE( "LogBlock", "[logging]" )
{
    CIE_TEST_CASE_INIT( "LogBlock" )

    std::filesystem::path loggerTestDir = getOutputPath();

    Logger logger(loggerTestDir / "LogBlock_test.log");

    for (Size i=0; i<10; ++i)
        if ( i % 2 == 0 )
        {
            LogBlock block( "iteration " + std::to_string(i), logger );
            block.log( "whatever" );
        }
}


} // namespace cie::utils