// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/logging/inc/Loggee.hpp"
#include "packages/logging/inc/Logger.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <iostream>
#include <string>
#include <vector>


namespace cie::utils {


std::filesystem::path loggeeTestDir = getOutputPath();


struct TestLoggee : public Loggee
{
    TestLoggee( Logger& logger ) : Loggee(logger,"TestLoggee") {}
};


CIE_TEST_CASE( "Loggee", "[logging]" )
{
    CIE_TEST_CASE_INIT( "Loggee" )

    // Create Logger
    Logger logger(loggeeTestDir / "loggee_test.txt");

    {
        // Create Loggee
        TestLoggee loggee( logger );
        size_t timerID;

        CIE_TEST_REQUIRE_NOTHROW( timerID = loggee.tic() );

        CIE_TEST_CHECK_NOTHROW( loggee.log( "test1" ) );
        CIE_TEST_CHECK_NOTHROW( loggee.log( "test2", LOG_TYPE_MESSAGE ) );
        CIE_TEST_CHECK_NOTHROW( loggee.log( "test3", LOG_TYPE_WARNING ) );
        CIE_TEST_CHECK_THROWS( loggee.log( "test4", LOG_TYPE_ERROR ) );

        CIE_TEST_CHECK_NOTHROW( loggee.logID( "test5", 0 ) );
        CIE_TEST_CHECK_NOTHROW( loggee.logID( "test6", 1, LOG_TYPE_MESSAGE ) );
        CIE_TEST_CHECK_NOTHROW( loggee.logID( "test7", 2, LOG_TYPE_WARNING ) );
        CIE_TEST_CHECK_THROWS( loggee.logID( "test8", 3, LOG_TYPE_ERROR ) );

        CIE_TEST_CHECK_NOTHROW( loggee.toc(timerID, false) );
        CIE_TEST_CHECK_NOTHROW( loggee.toc( "test9", timerID, false) );
        CIE_TEST_CHECK_NOTHROW( loggee.toc(timerID) );

        CIE_TEST_CHECK_NOTHROW( loggee << "test10" );
        CIE_TEST_CHECK_NOTHROW( loggee << 11 );
        CIE_TEST_CHECK_NOTHROW( loggee << std::vector<int>{ 12, 13 } );

        CIE_TEST_CHECK_NOTHROW( loggee.separate() );
    }
}


}