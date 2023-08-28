#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/Plot2.hpp"

// --- STL Includes ---
#include <deque>
#include <vector>
#include <math.h>
#include <cmath>

// --- STL Includes ---
#include <chrono>
#include <thread>


namespace cie::gl {


CIE_TEST_CASE( "Plot2", "[plot]" )
{
    CIE_TEST_CASE_INIT( "Plot2" )

    using XContainer = std::deque<Double>;
    using YContainer = std::vector<Double>;
    XContainer X;
    YContainer Y;

    Size n    = 100;
    double dx = 2.0 * M_PI / double(n);

    for ( double x=0.0; x<=2.0*M_PI; x+=dx )
    {
        X.push_back( x );
        Y.push_back( std::cos(x) );
    }

    Plot2<XContainer,YContainer> plot( X, Y );
    plot.update();
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    //plot2( X, Y );
}


} // namespace cie::gl