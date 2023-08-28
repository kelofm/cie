// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/logging/inc/Logger.hpp"
#include "cmake_variables.hpp"


namespace cie::utils {


CIE_TEST_CASE( "Logger", "[logger]" )
{
    CIE_TEST_CASE_INIT( "Logger" )

    std::filesystem::path loggerTestDir = getOutputPath();

    {
        // Construct logger
        Logger logger(loggerTestDir / "testLog.txt");

        CIE_TEST_CHECK_NOTHROW( logger.log( "test1" ) );
        CIE_TEST_CHECK_NOTHROW( logger.warn( "test2" ) );
        CIE_TEST_CHECK_THROWS( logger.error( "test3" ) );
        CIE_TEST_CHECK_NOTHROW( logger.separate( ) );
        CIE_TEST_CHECK_NOTHROW( logger.increaseIndent() );
        CIE_TEST_CHECK_NOTHROW( logger.log( "test4" ) );
        CIE_TEST_CHECK_NOTHROW( logger.increaseIndent() );
        CIE_TEST_CHECK_NOTHROW( logger.increaseIndent() );
        CIE_TEST_CHECK_NOTHROW( logger.log( "test5" ) );
        CIE_TEST_CHECK_NOTHROW( logger.decreaseIndent() );
        CIE_TEST_CHECK_NOTHROW( logger.log( "test6" ) );
        size_t timerID = 0;
        CIE_TEST_CHECK_NOTHROW( timerID = logger.startTimer() );
        CIE_TEST_CHECK_NOTHROW( logger.noIndent() );
        CIE_TEST_CHECK_NOTHROW( logger.logElapsed( "test7", timerID ) );
    }
}


}