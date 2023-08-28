// --- Internal Includes ---
#include "packages/utilities/inc/errortypes.hpp"
#include "cmake_variables.hpp"


namespace cie::linalg {

MatrixError::MatrixError( const std::string& r_location,
                          const std::string& r_message ) :
    Exception( r_location, "Matrix error: " + r_message )
{
}

} // namespace cie::linalg