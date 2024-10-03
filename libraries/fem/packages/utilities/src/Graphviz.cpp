// --- FEM Includes ---
#include "packages/utilities/inc/Graphviz.hpp"

namespace cie::io {


Graphviz::Output::Output(Ref<std::ostream> rStream, Settings settings)
    : _pStream(&rStream),
      _settings(settings)
{
}


Graphviz::Output::Output()
    : Output(std::cout, Settings())
{
}





} // namespace cie::io
