// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/control/inc/utilities.hpp"


namespace cie::gl {


void transformScreenCoordinates( double& x, double& y, AbsWindow::WeakPointer wp_window )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( auto p_window = wp_window.lock() )
    {
        const auto& r_windowSize = p_window->getSize();

        CIE_DIVISION_BY_ZERO_CHECK( r_windowSize.first != 0 )
        CIE_DIVISION_BY_ZERO_CHECK( r_windowSize.second != 0 )

        x /= r_windowSize.first;
        y = ( r_windowSize.second - y ) / r_windowSize.second;
    }
    else
        CIE_THROW( Exception, "Pointer to window expired!" )

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl