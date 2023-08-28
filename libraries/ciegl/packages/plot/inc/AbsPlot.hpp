#ifndef CIE_GL_PLOT_ABS_PLOT_HPP
#define CIE_GL_PLOT_ABS_PLOT_HPP

// --- Internal Includes ---
#include "packages/context/inc/AbsContext.hpp"
#include "packages/context/inc/AbsWindow.hpp"
#include "packages/buffer/inc/AttributeContainer.hpp"

// --- STL Includes ---
#include <vector>
#include <deque>
#include <memory>


namespace cie::gl {


class AbsPlot
{
public:
    AbsPlot( AbsWindow::SharedPointer p_window );
    AbsPlot();

    virtual void update();
    virtual void fit( bool keepAspectRatio = false ) = 0;

    void show();

    const AttributeContainer& attributes() const;

protected:
    AbsWindow::SharedPointer _p_window;
    AttributeContainer::SharedPointer _p_attributes;
};


} // namespace cie::gl


#endif