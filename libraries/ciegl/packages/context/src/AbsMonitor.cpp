// --- Internal Includes ---
#include "packages/context/inc/AbsMonitor.hpp"


namespace cie::gl {


AbsMonitor::AbsMonitor( Size id,
                        Size width,
                        Size height,
                        bool isPrimary ) :
    utils::IDObject<Size>(id),
    _resolution( width, height ),
    _isPrimary( isPrimary )
{
}


AbsMonitor::~AbsMonitor()
{
}


const std::pair<Size,Size>& AbsMonitor::resolution() const
{
    return _resolution;
}


bool AbsMonitor::isPrimary() const
{
    return _isPrimary;
}


} // namespace cie::gl