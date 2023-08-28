// --- Internal Includes ---
#include "packages/control/inc/MappedCameraControls.hpp"


namespace cie::gl {


MappedCameraControls::MappedCameraControls() :
    AbsCameraControls(),
    MappedControls()
{
}


void MappedCameraControls::onSubjectChange()
{
    for ( auto key : this->_activeKeys )
        this->callControl( key );
}


} // namespace cie::gl